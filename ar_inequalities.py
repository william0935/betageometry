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
"""

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from relations import Point, RelationNode

# Coefficients below this are treated as zero when reading a solution back.
_EPS = 1e-9


def nonneg_combination(nonneg_rows: Sequence[Sequence[float]],
                       free_rows: Sequence[Sequence[float]],
                       target: Sequence[float],
                       tol: float = 1e-9) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Express `target` as a combination of the given rows, or return None.

    The coefficients on `nonneg_rows` are constrained to be >= 0; those on `free_rows`
    may take either sign. Returns `(lambdas, mus)` or None if no such combination
    exists.

    This is a linear feasibility problem -- find z >= 0 with M z = target, where the
    free rows appear as a difference of two non-negative columns -- solved with a
    phase-1 simplex. Bland's rule is used for pivot selection: it is slower than a
    steepest-edge rule but cannot cycle, and these systems are tiny.
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

    solution = _phase_one_simplex(matrix, rhs, tol=tol)
    if solution is None:
        return None

    lambdas = solution[:k_nonneg]
    mus = solution[k_nonneg:k_nonneg + k_free] - solution[k_nonneg + k_free:]
    return lambdas, mus


def _phase_one_simplex(matrix: np.ndarray, rhs: np.ndarray,
                       tol: float = 1e-9, max_iter: int = 10000) -> Optional[np.ndarray]:
    """Find z >= 0 with `matrix @ z == rhs`, or None if none exists."""
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

    for _ in range(max_iter):
        # Bland's rule: the lowest-indexed column with a negative reduced cost.
        entering = -1
        for j in range(n + m):
            if cost[j] < -tol:
                entering = j
                break
        if entering == -1:
            break  # optimal

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
            break  # unbounded along this column, which cannot happen in phase 1

        pivot = tableau[leaving, entering]
        tableau[leaving, :] /= pivot
        for i in range(m):
            if i != leaving and abs(tableau[i, entering]) > tol:
                tableau[i, :] -= tableau[i, entering] * tableau[leaving, :]
        if abs(cost[entering]) > tol:
            cost -= cost[entering] * tableau[leaving, :]
        basis[leaving] = entering

    # The objective value is -cost[-1]; feasible exactly when every artificial is zero.
    if -cost[n + m] > tol * max(1.0, float(np.abs(b).sum())):
        return None

    solution = np.zeros(n + m)
    for i, var in enumerate(basis):
        solution[var] = tableau[i, -1]
    if np.any(solution[n:] > tol):
        return None  # an artificial stayed in the basis at a non-zero level
    return solution[:n]


class InequalityTable:
    """Linear inequalities over a set of named columns.

    Every stored row means ``coefficients . x >= 0``, or ``> 0`` when `strict`.
    Equalities (``== 0``) are stored separately, since their coefficients may be scaled
    in either direction when combining.
    """

    def __init__(self, header: Optional[List[Any]] = None):
        self.header: List[Any] = list(header) if header else []
        self.col_id: Dict[Any, int] = {c: i for i, c in enumerate(self.header)}
        # Parallel lists: coefficients, strictness, and the relation that supplied it.
        self.ineq_rows: List[List[float]] = []
        self.ineq_strict: List[bool] = []
        self.ineq_source: List[RelationNode] = []
        self.eq_rows: List[List[float]] = []
        self.eq_source: List[RelationNode] = []

    def add_col(self, name: Any) -> int:
        if name not in self.col_id:
            self.col_id[name] = len(self.header)
            self.header.append(name)
        return self.col_id[name]

    def table_length(self) -> int:
        return len(self.header)

    def _widen(self, row: Sequence[float]) -> List[float]:
        row = list(row)
        return row + [0.0] * (len(self.header) - len(row))

    def add_inequality(self, row: Sequence[float], relation: RelationNode,
                       strict: bool = True):
        """Record ``row . x >= 0`` (or ``> 0``).

        Unlike `ar.Table`, a row that is linearly dependent on existing rows is still
        kept: `a >= b` alongside `b >= a` is dependent but carries real information.
        """
        self.ineq_rows.append(self._widen(row))
        self.ineq_strict.append(strict)
        self.ineq_source.append(relation)

    def add_equality(self, row: Sequence[float], relation: RelationNode):
        """Record ``row . x == 0``, usable with a coefficient of either sign."""
        self.eq_rows.append(self._widen(row))
        self.eq_source.append(relation)

    def _padded(self, rows: List[List[float]]) -> List[List[float]]:
        width = len(self.header)
        return [r + [0.0] * (width - len(r)) for r in rows]

    def implies(self, row: Sequence[float], strict: bool = True,
                tol: float = 1e-9) -> Tuple[bool, frozenset]:
        """Is ``row . x >= 0`` (or ``> 0``) a consequence of what is stored?

        Returns `(proved, sources)` to match `ar.Table.is_spanned`, so callers and the
        proof trace treat both tables the same way.
        """
        if not self.ineq_rows and not self.eq_rows:
            return False, frozenset()

        target = self._widen(row)
        nonneg = self._padded(self.ineq_rows)
        free = self._padded(self.eq_rows)

        result = nonneg_combination(nonneg, free, target, tol=tol)
        if result is None:
            return False, frozenset()

        lambdas, mus = result
        used_ineq = [i for i, v in enumerate(lambdas) if v > _EPS]

        # A strict conclusion needs at least one strict premise pulling its weight. A
        # sum of non-strict inequalities can only ever prove a non-strict one.
        if strict and not any(self.ineq_strict[i] for i in used_ineq):
            return False, frozenset()

        sources = {self.ineq_source[i] for i in used_ineq}
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
