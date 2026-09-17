"""Algebraic reasoning over inequalities.

Deliberately separate from `ar.py` rather than a subclass of `Table`, for two reasons.

**Only addition is allowed.** `Table` asks whether a row lies in the *span* of its rows,
which permits any real coefficient. Scaling an inequality by a negative number reverses
it, so a span test concludes `CD >= AB` from `AB >= CD`. Here a row is implied only by a
combination with **non-negative** coefficients on the inequality rows (equalities may
still be scaled freely, in either direction). That is Farkas' lemma, and it is the whole
difference between the two modules.

**The column space is different.** `AngleTable`'s columns are line directions modulo pi,
which is correct for equalities and meaningless for inequalities -- mod pi there is no
order, since adding 180 degrees flips any comparison. So the angle inequality table has
one column per *angle magnitude* (a vertex and the two points it opens towards) rather
than per direction.

Rows are kept even when linearly dependent on rows already present. `Table` drops those
as redundant, which is right for equalities and lossy here: `a >= b` and `b >= a` are
linearly dependent but together prove `a == b`.

Strict (`>`) and non-strict (`>=`) rows are held in separate tables. A strict conclusion
needs a combination that puts real weight on a strict premise -- a sum of non-strict rows
can only ever prove a non-strict claim -- but the two are not otherwise independent, and
querying the strict table alone would be wrong: `a > b` together with `b >= c` is a
perfectly good strict proof of `a > c`, and it draws a row from each. So the split is not
two separate searches but one search that is *told* to lean on the strict table, and the
distinction is enforced while choosing the certificate rather than checked afterwards.
Checking afterwards loses real derivations: with `a >= b`, `b >= c` and `a > c` all known,
a solver free to return any certificate may answer `a >= b >= c`, which uses no strict row,
and `a > c` would be reported unprovable although it was a premise.
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from relations import Point, RelationNode

# Coefficients below this are treated as zero when reading a solution back.
_EPS = 1e-9


def nonneg_combination(nonneg_rows: Sequence[Sequence[float]],
                       free_rows: Sequence[Sequence[float]],
                       target: Sequence[float],
                       tol: float = 1e-9,
                       weights: Optional[Sequence[float]] = None
                       ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Express `target` as a combination of the given rows, or return None.

    The coefficients on `nonneg_rows` are constrained to be >= 0; those on `free_rows`
    may take either sign. Returns `(lambdas, mus)` or None if no such combination
    exists.

    This is a linear feasibility problem -- find z >= 0 with M z = target, where the
    free rows appear as a difference of two non-negative columns -- solved with a
    phase-1 simplex. Bland's rule is used for pivot selection: it is slower than a
    steepest-edge rule but cannot cycle, and these systems are tiny.

    With `weights`, one per non-negative row, the search does not stop at the first
    combination it finds: a phase-2 pass then maximises `weights . lambdas`, so among all
    the certificates that exist the one returned carries as much weight as possible on
    the rows the caller cares about. `InequalityTable.implies` uses this to ask for a
    certificate that actually rests on a strict premise.
    """
    n_cols = len(target)
    k_nonneg = len(nonneg_rows)
    k_free = len(free_rows)
    if k_nonneg == 0 and k_free == 0:
        return (np.zeros(0), np.zeros(0)) if not any(target) else None

    # Columns: [ nonneg rows | free rows | -free rows ]
    columns = []
    for row in nonneg_rows:
        columns.append(np.asarray(row, dtype=float))
    for row in free_rows:
        columns.append(np.asarray(row, dtype=float))
    for row in free_rows:
        columns.append(-np.asarray(row, dtype=float))

    matrix = np.column_stack(columns) if columns else np.zeros((n_cols, 0))
    rhs = np.asarray(target, dtype=float)

    objective = None
    if weights is not None:
        objective = np.zeros(matrix.shape[1])
        objective[:k_nonneg] = np.asarray(weights, dtype=float)

    solution = _simplex(matrix, rhs, objective=objective, tol=tol)
    if solution is None:
        return None

    lambdas = solution[:k_nonneg]
    mus = solution[k_nonneg:k_nonneg + k_free] - solution[k_nonneg + k_free:]
    return lambdas, mus


def _pivot(tableau: np.ndarray, basis: List[int], row: int, col: int):
    """Make `col` basic in `row`, by Gauss-Jordan on the tableau."""
    tableau[row, :] /= tableau[row, col]
    for i in range(tableau.shape[0]):
        if i != row and tableau[i, col] != 0.0:
            tableau[i, :] -= tableau[i, col] * tableau[row, :]
    basis[row] = col


def _pivot_to_optimal(tableau: np.ndarray, basis: List[int], cost: np.ndarray,
                      n_cols: int, tol: float, max_iter: int) -> int:
    """Minimise `cost` over the first `n_cols` columns, in place.

    Returns -1 once no reduced cost is negative, or the index of a column along which
    the objective falls without bound.
    """
    m = tableau.shape[0]
    for _ in range(max_iter):
        # Bland's rule: the lowest-indexed column with a negative reduced cost.
        entering = -1
        for j in range(n_cols):
            if cost[j] < -tol:
                entering = j
                break
        if entering == -1:
            return -1  # optimal

        column = tableau[:, entering]
        # Ratio test, breaking ties on the lowest basis index (Bland).
        leaving, best_ratio = -1, None
        for i in range(m):
            if column[i] > tol:
                ratio = tableau[i, -1] / column[i]
                if best_ratio is None or ratio < best_ratio - tol or \
                   (abs(ratio - best_ratio) <= tol and basis[i] < basis[leaving]):
                    leaving, best_ratio = i, ratio
        if leaving == -1:
            return entering  # unbounded along this column

        coeff = cost[entering]
        _pivot(tableau, basis, leaving, entering)
        if abs(coeff) > tol:
            cost -= coeff * tableau[leaving, :]
    return -1


def _phase_one(matrix: np.ndarray, rhs: np.ndarray, tol: float, max_iter: int):
    """Find a feasible basis for z >= 0, `matrix @ z == rhs`.

    Returns `(tableau, basis, n)` with the artificial columns already removed, or None
    if the system is infeasible.
    """
    m, n = matrix.shape
    a = matrix.astype(float).copy()
    b = rhs.astype(float).copy()

    # Artificial variables need a non-negative right-hand side.
    for i in range(m):
        if b[i] < 0:
            a[i, :] *= -1
            b[i] *= -1

    # Tableau with one artificial per row; the artificials form the initial basis.
    tableau = np.hstack([a, np.eye(m), b.reshape(-1, 1)])
    basis = list(range(n, n + m))

    # Minimise the sum of artificials. The objective row starts as minus the sum of the
    # constraint rows, because every artificial is basic at the outset.
    cost = np.zeros(n + m + 1)
    cost[:n] = -a.sum(axis=0)
    cost[n + m] = -b.sum()

    _pivot_to_optimal(tableau, basis, cost, n + m, tol, max_iter)

    # The objective value is -cost[-1]; feasible exactly when every artificial is zero.
    if -cost[n + m] > tol * max(1.0, float(np.abs(b).sum())):
        return None

    # An artificial may still be basic at a zero level. Pivot each one out on any
    # structural column with a non-zero entry in its row; a row with none is redundant
    # and is dropped. Leaving one basic would let phase 2 raise it off zero, which
    # silently breaks feasibility -- the ratio test only guards rows whose entry in the
    # entering column is positive.
    redundant = []
    for i in range(m):
        if basis[i] < n:
            continue
        col = -1
        for j in range(n):
            if abs(tableau[i, j]) > tol:
                col = j
                break
        if col == -1:
            redundant.append(i)
        else:
            _pivot(tableau, basis, i, col)
    if redundant:
        drop = set(redundant)
        keep = [i for i in range(m) if i not in drop]
        tableau = tableau[keep, :]
        basis = [basis[i] for i in keep]

    # Drop the artificial columns; nothing downstream may reintroduce them.
    tableau = np.hstack([tableau[:, :n], tableau[:, -1:]])
    return tableau, basis, n


def _simplex(matrix: np.ndarray, rhs: np.ndarray,
             objective: Optional[np.ndarray] = None,
             tol: float = 1e-9, max_iter: int = 10000) -> Optional[np.ndarray]:
    """Find z >= 0 with `matrix @ z == rhs`, or None if none exists.

    With `objective`, the solution returned maximises `objective . z` rather than being
    whichever feasible point phase 1 happened to stop at.
    """
    found = _phase_one(matrix, rhs, tol=tol, max_iter=max_iter)
    if found is None:
        return None
    tableau, basis, n = found

    extra = None
    if objective is not None:
        extra = _phase_two(tableau, basis, n, np.asarray(objective, dtype=float),
                           tol=tol, max_iter=max_iter)

    solution = np.zeros(n)
    for i, var in enumerate(basis):
        if var < n:
            solution[var] = tableau[i, -1]
    if extra is not None:
        col, value = extra
        solution[col] += value
    return solution


def _phase_two(tableau: np.ndarray, basis: List[int], n: int, objective: np.ndarray,
               tol: float, max_iter: int) -> Optional[Tuple[int, float]]:
    """Maximise `objective . z` from a feasible basis, in place.

    Returns None normally, or `(column, value)` for a non-basic variable that has to be
    added to the basic solution when the objective turned out to be unbounded.
    """
    # The pivot loop minimises, so it is handed the negated objective.
    d = -objective
    cost = np.zeros(n + 1)
    cost[:n] = d
    for i, var in enumerate(basis):
        if d[var] != 0.0:
            cost -= d[var] * tableau[i, :]

    unbounded_col = _pivot_to_optimal(tableau, basis, cost, n, tol, max_iter)
    if unbounded_col < 0:
        return None

    # The objective grows without bound along this column, so any positive step gives a
    # certificate worth more than the current one. A single unit is enough: the caller
    # only needs to know the weight it asked about can be made positive. Every entry of
    # the column is <= 0, so the basic variables only increase and stay feasible.
    step = 1.0
    tableau[:, -1] -= step * tableau[:, unbounded_col]
    return unbounded_col, step


class InequalityTable:
    """Linear inequalities over a set of named columns.

    Every stored row means ``coefficients . x >= 0``, or ``> 0`` when `strict`. Strict and
    non-strict rows are kept in separate tables so that a strict query can require the
    certificate to draw on the strict one; see the module docstring. Equalities (``== 0``)
    are separate again, since their coefficients may be scaled in either direction when
    combining.
    """

    def __init__(self, header: Optional[List[Any]] = None):
        self.header: List[Any] = list(header) if header else []
        self.col_id: Dict[Any, int] = {c: i for i, c in enumerate(self.header)}
        # Parallel lists of coefficients and the relation that supplied each row.
        self.strict_rows: List[List[float]] = []
        self.strict_source: List[RelationNode] = []
        self.weak_rows: List[List[float]] = []
        self.weak_source: List[RelationNode] = []
        self.eq_rows: List[List[float]] = []
        self.eq_source: List[RelationNode] = []

    def add_col(self, name: Any) -> int:
        if name not in self.col_id:
            self.col_id[name] = len(self.header)
            self.header.append(name)
        return self.col_id[name]

    def table_length(self) -> int:
        return len(self.header)

    def clear_rows(self):
        """Drop every row, keeping the columns.

        `DDWithAR` reseeds the tables after each deduction pass; the header is left alone
        so column ids stay stable across passes.
        """
        self.strict_rows, self.strict_source = [], []
        self.weak_rows, self.weak_source = [], []
        self.eq_rows, self.eq_source = [], []

    def _widen(self, row: Sequence[float]) -> List[float]:
        row = list(row)
        return row + [0.0] * (len(self.header) - len(row))

    def _padded(self, rows: List[List[float]]) -> List[List[float]]:
        """Rows stored before the header grew still have to reach its full width."""
        return [self._widen(r) for r in rows]

    def add_inequality(self, row: Sequence[float], relation: RelationNode,
                       strict: bool = True):
        """Record ``row . x >= 0`` (or ``> 0``), in whichever table matches.

        Unlike `ar.Table`, a row that is linearly dependent on existing rows is still
        kept: `a >= b` alongside `b >= a` is dependent but carries real information.
        """
        if strict:
            self.strict_rows.append(self._widen(row))
            self.strict_source.append(relation)
        else:
            self.weak_rows.append(self._widen(row))
            self.weak_source.append(relation)

    def add_equality(self, row: Sequence[float], relation: RelationNode):
        """Record ``row . x == 0``, usable with a coefficient of either sign."""
        self.eq_rows.append(self._widen(row))
        self.eq_source.append(relation)

    def implies(self, row: Sequence[float], strict: bool = True,
                tol: float = 1e-9) -> Tuple[bool, frozenset]:
        """Is ``row . x >= 0`` (or ``> 0``) a consequence of what is stored?

        Returns `(proved, sources)` to match `ar.Table.is_spanned`, so callers and the
        proof trace treat both tables the same way.
        """
        if not (self.strict_rows or self.weak_rows or self.eq_rows):
            return False, frozenset()
        if strict and not self.strict_rows:
            return False, frozenset()  # nothing that could make the conclusion strict

        target = self._widen(row)
        # A strict row serves as a non-strict one too -- `a > b` gives `a >= b` -- so both
        # tables go to the search. It is the weighting below, not the choice of rows, that
        # separates a strict conclusion from a non-strict one.
        n_strict = len(self.strict_rows)
        nonneg = self._padded(self.strict_rows) + self._padded(self.weak_rows)
        free = self._padded(self.eq_rows)

        # For a strict goal, ask for the certificate leaning hardest on the strict table;
        # if even that puts no weight there, no certificate does. See the module docstring
        # for why this cannot be decided after the certificate is chosen.
        weights = [1.0] * n_strict + [0.0] * len(self.weak_rows) if strict else None

        result = nonneg_combination(nonneg, free, target, tol=tol, weights=weights)
        if result is None:
            return False, frozenset()

        lambdas, mus = result
        if strict and float(np.sum(lambdas[:n_strict])) <= _EPS:
            return False, frozenset()

        sources = set()
        for i, value in enumerate(lambdas):
            if value > _EPS:
                sources.add(self.strict_source[i] if i < n_strict
                            else self.weak_source[i - n_strict])
        sources |= {self.eq_source[j] for j, v in enumerate(mus) if abs(v) > _EPS}
        return True, frozenset(s for s in sources if s is not None)


def segment_key(p1: Point, p2: Point) -> frozenset:
    """Column key for the length |p1p2|."""
    return frozenset({p1, p2})


def angle_key(p1: Point, vertex: Point, p2: Point) -> Tuple[Point, frozenset]:
    """Column key for the magnitude of angle p1-vertex-p2.

    Keyed on the vertex plus the unordered pair it opens towards, mirroring how the rest
    of the codebase keys segments and circles. Note this is a magnitude in [0, 180], not
    a direction modulo pi -- ordering only makes sense for the former.
    """
    return (vertex, frozenset({p1, p2}))


class SegmentInequalityTable(InequalityTable):
    """Inequalities between segment lengths."""

    def ensure_segment(self, p1: Point, p2: Point) -> int:
        return self.add_col(segment_key(p1, p2))

    def row_for_difference(self, a1: Point, a2: Point, b1: Point, b2: Point) -> List[float]:
        """A row meaning ``|a1a2| - |b1b2|``."""
        self.ensure_segment(a1, a2)
        self.ensure_segment(b1, b2)
        row = [0.0] * len(self.header)
        row[self.col_id[segment_key(a1, a2)]] += 1.0
        row[self.col_id[segment_key(b1, b2)]] -= 1.0
        return row


class AngleInequalityTable(InequalityTable):
    """Inequalities between angle magnitudes."""

    def ensure_angle(self, p1: Point, vertex: Point, p2: Point) -> int:
        return self.add_col(angle_key(p1, vertex, p2))

    def row_for_difference(self, a1: Point, av: Point, a2: Point,
                           b1: Point, bv: Point, b2: Point) -> List[float]:
        """A row meaning ``angle(a1 av a2) - angle(b1 bv b2)``."""
        self.ensure_angle(a1, av, a2)
        self.ensure_angle(b1, bv, b2)
        row = [0.0] * len(self.header)
        row[self.col_id[angle_key(a1, av, a2)]] += 1.0
        row[self.col_id[angle_key(b1, bv, b2)]] -= 1.0
        return row
