import math
from relations import *
import numpy as np
from typing import List, Tuple, Dict, Any, Optional

# checks sameclock
def is_sameclock(p1: Point, p2: Point, p3: Point, p4: Point, p5: Point, p6: Point) -> bool:
    return ((p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x)) * \
            ((p5.x - p4.x) * (p6.y - p4.y) - (p5.y - p4.y) * (p6.x - p4.x)) > 0


# Special column key: represents a fixed +90 degrees constant in column 0
CONST_90 = "90_degrees"

# Entries below this are treated as exact zeros during elimination. Rows are built from
# small integers, so anything at this magnitude is floating-point noise.
_PIVOT_EPS = 1e-14


class Table:
    """A set of linear constraints over segment/angle columns.

    Membership ("is this row in the span of the rows I already have?") is the single
    hottest operation in the solver, so the table keeps an incremental row-echelon
    basis rather than re-factorising from scratch on every query. A query is one
    O(rows x cols) elimination sweep; the elimination coefficients double as the
    dependency certificate, so no separate least-squares solve is needed.

    Every row that survives `add_row` is by construction linearly independent of the
    rows before it, so there is exactly one basis vector per stored row.
    """

    def __init__(self, header: List[Any]):
        self.header = list(header)
        self.col_id: Dict[Any, int] = {col: i for i, col in enumerate(self.header)}
        # Kept as plain lists purely so the HTML dump can render the original rows.
        self.rows: Dict[RelationNode, List[List[int]]] = {}
        self.relations: set[RelationNode] = set()  # Track all relations in order

        # Echelon basis, stored in fixed-capacity arrays that double on demand so that
        # adding a column or a row is amortised O(1) rather than a full reallocation.
        self._n = 0                                  # number of basis vectors / stored rows
        self._basis = np.zeros((8, 8))               # _basis[i][_pivots[i]] == 1
        self._coeff = np.zeros((8, 8))               # _basis[i] == sum_k _coeff[i][k] * row_k
        self._pivots: List[int] = []
        self._row_relation: List[RelationNode] = []
        # Memo for is_spanned. The rules issue the same queries over and over -- the
        # triangle rules alone ask about every shared segment of every candidate pair,
        # on every pass.
        #
        # `_span_cache` is never invalidated: `add_row` only ever stores rows that are
        # outside the current span, so the span grows monotonically and a row that is
        # spanned stays spanned. The stored rows also stay linearly independent, so the
        # combination that expresses a spanned row -- the dependency certificate -- is
        # unique and does not change when later rows arrive.
        #
        # A miss, on the other hand, can turn into a hit as soon as a row is added, so
        # `_miss_cache` is dropped whenever the table changes.
        self._span_cache: Dict[tuple, Tuple[bool, frozenset]] = {}
        self._miss_cache: set = set()

    def _ensure_capacity(self, n_rows: int, n_cols: int):
        cap_r, cap_c = self._basis.shape
        new_r, new_c = cap_r, cap_c
        while new_r < n_rows:
            new_r *= 2
        while new_c < n_cols:
            new_c *= 2
        if new_r == cap_r and new_c == cap_c:
            return
        basis = np.zeros((new_r, new_c))
        basis[:cap_r, :cap_c] = self._basis
        self._basis = basis
        coeff = np.zeros((new_r, new_r))
        coeff[:cap_r, :cap_r] = self._coeff
        self._coeff = coeff

    def _eliminate(self, row: List[Any]) -> Tuple[np.ndarray, np.ndarray]:
        """Reduce `row` against the basis.

        Returns (residual, comb) with  row == residual + sum_k comb[k] * row_k.
        """
        n = self._n
        work = np.zeros(self._basis.shape[1])
        work[:len(row)] = row
        comb = np.zeros(self._basis.shape[0])
        for i in range(n):
            f = work[self._pivots[i]]
            if abs(f) > _PIVOT_EPS:
                work -= f * self._basis[i]
                comb[:n] += f * self._coeff[i, :n]
        return work, comb

    def add_col(self, col_name: Any):
        if col_name in self.col_id:
            return
        self.header.append(col_name)
        self.col_id[col_name] = len(self.header) - 1
        # Extend the display rows; the basis is zero-padded by construction, so it
        # needs no work here beyond having room for the new column.
        for row_list in self.rows.values():
            for row in row_list:
                row.append(0)
        self._ensure_capacity(self._n, len(self.header))
        self._miss_cache.clear()

    def is_spanned(self, row: List[Any], tol: float = 1e-10) -> Tuple[bool, frozenset]:
        # The zero row is the empty constraint ("0 = 0"), which every set of rows
        # implies, including none at all. It needs no premises to justify it.
        if not any(row):
            return True, frozenset()
        if self._n == 0:
            return False, frozenset()

        key = tuple(row)
        cached = self._span_cache.get(key)
        if cached is not None:
            return cached
        if key in self._miss_cache:
            return False, frozenset()

        work, comb = self._eliminate(row)
        if np.abs(work).max() > tol:
            self._miss_cache.add(key)
            return False, frozenset()

        # frozenset so the memoised value cannot be mutated through a caller
        result = (True, frozenset(self._row_relation[k]
                                  for k in np.nonzero(np.abs(comb[:self._n]) > tol)[0]))
        self._span_cache[key] = result
        return result

    def add_row(self, row: List[Any], relation: RelationNode, tol: float = 1e-10):
        if len(row) != len(self.header):
            raise ValueError("Row length does not match header length.")

        n = self._n
        self._ensure_capacity(n + 1, len(self.header))
        work, comb = self._eliminate(row)

        nonzero = np.abs(work) > tol
        if not nonzero.any():
            return  # already implied by the existing rows

        self._miss_cache.clear()

        # row_n == work + sum_k comb[k] * row_k, so normalising work by its leading
        # entry gives a basis vector expressed over the stored rows as (e_n - comb)/scale.
        pivot = int(np.argmax(nonzero))
        scale = work[pivot]
        self._basis[n] = work / scale
        self._coeff[n] = -comb / scale
        self._coeff[n, n] += 1.0 / scale
        self._pivots.append(pivot)
        self._row_relation.append(relation)
        self._n = n + 1

        if relation not in self.rows:
            self.rows[relation] = []
        self.rows[relation].append(list(row))
        self.relations.add(relation)

    def table_length(self) -> int:
        return len(self.header)

class AngleTable(Table):
    def __init__(self, header: List[Any]):
        base_header = [col for col in header if col != CONST_90]
        ordered_header = [CONST_90] + base_header
        super().__init__(ordered_header)
    
    def add_eqangle(self, eqangle: EqualAngle):
        p1, p2, p3, p4, p5, p6 = eqangle.points
        seg1_1 = frozenset({p1, p2})
        seg1_2 = frozenset({p2, p3})
        seg2_1 = frozenset({p4, p5})
        seg2_2 = frozenset({p5, p6})
        
        coefficients = {}
        coefficients[seg1_1] = coefficients.get(seg1_1, 0) + 1
        coefficients[seg1_2] = coefficients.get(seg1_2, 0) - 1
        coefficients[seg2_1] = coefficients.get(seg2_1, 0) - 1
        coefficients[seg2_2] = coefficients.get(seg2_2, 0) + 1
        for seg in coefficients:
            if seg not in self.col_id:
                self.add_col(seg)
        
        row = [0] * len(self.header)
        for seg, coeff in coefficients.items():
            if coeff != 0:
                row[self.col_id[seg]] = coeff

        self.add_row(row, eqangle)
    
    def add_parallel(self, parallel: Parallel):
        p1, p2, p3, p4 = parallel.points
        seg1 = frozenset({p1, p2})
        seg2 = frozenset({p3, p4})

        coefficients = {}
        coefficients[seg1] = coefficients.get(seg1, 0) + 1
        coefficients[seg2] = coefficients.get(seg2, 0) - 1
        for seg in coefficients:
            if seg not in self.col_id:
                self.add_col(seg)
        
        row = [0] * len(self.header)
        for seg, coeff in coefficients.items():
            if coeff != 0:
                row[self.col_id[seg]] = coeff

        self.add_row(row, parallel)

    def add_collinear(self, collinear : Collinear):
        p1, p2, p3 = collinear.points
        s1 = frozenset({p1, p2})
        s2 = frozenset({p2, p3})
        s3 = frozenset({p1, p3})
        
        def helper_add_collinear(s1, s2):
            coefficients = {}
            coefficients[s1] = coefficients.get(s1, 0) + 1
            coefficients[s2] = coefficients.get(s2, 0) - 1
        
            for seg in coefficients:
                if seg not in self.col_id:
                    self.add_col(seg)
            
            row = [0] * len(self.header)
            for seg, coeff in coefficients.items():
                if coeff != 0:
                    row[self.col_id[seg]] = coeff

            self.add_row(row, collinear)

        helper_add_collinear(s1, s2)
        helper_add_collinear(s2, s3)

    def add_perpendicular(self, perpendicular: Perpendicular):
        p1, p2, p3, p4 = perpendicular.points
        seg1 = frozenset({p1, p2})
        seg2 = frozenset({p3, p4})

        # The row encodes dir(seg1) - dir(seg2) - 90 = 0, so seg1 must be the line whose
        # direction is the larger of the two. The tie-break only has to be *consistent*
        # for a given pair of segments -- picking the opposite orientation for the same
        # pair anywhere else would let the solver derive CONST_90 == 0 and collapse the
        # table -- which is why direction_mod_pi pins the pi end of the range to 0.
        angle1 = direction_mod_pi(p1, p2)
        angle2 = direction_mod_pi(p3, p4)

        if angle1 < angle2:
            seg1, seg2 = seg2, seg1

        coefficients: Dict[Any, int] = {}
        coefficients[seg1] = coefficients.get(seg1, 0) + 1
        coefficients[seg2] = coefficients.get(seg2, 0) - 1
        coefficients[CONST_90] = coefficients.get(CONST_90, 0) - 1

        for seg in (seg1, seg2):
            if seg not in self.col_id:
                self.add_col(seg)

        row = [0] * len(self.header)
        for key, coeff in coefficients.items():
            row[self.col_id[key]] = coeff

        self.add_row(row, perpendicular)

class RatioTable(Table):
    def __init__(self, header: List[frozenset[Point]]):
        super().__init__(header)
    
    def add_eqratio(self, eqratio: EqualRatio):
        p1, p2, p3, p4, p5, p6, p7, p8 = eqratio.points

        coefficients = {}
        s12 = frozenset({p1, p2})
        s34 = frozenset({p3, p4})
        s56 = frozenset({p5, p6})
        s78 = frozenset({p7, p8})
        coefficients[s12] = coefficients.get(s12, 0) + 1
        coefficients[s34] = coefficients.get(s34, 0) - 1
        coefficients[s56] = coefficients.get(s56, 0) - 1
        coefficients[s78] = coefficients.get(s78, 0) + 1

        for seg in coefficients:
            if seg not in self.col_id:
                self.add_col(seg)

        row = [0] * len(self.header)
        for seg, coeff in coefficients.items():
            if coeff != 0:
                row[self.col_id[seg]] = coeff

        self.add_row(row, eqratio)

    def add_cong(self, cong: Congruent):
        p1, p2, p3, p4 = cong.points

        coefficients = {}
        s12 = frozenset({p1, p2})
        s34 = frozenset({p3, p4})
        coefficients[s12] = coefficients.get(s12, 0) + 1
        coefficients[s34] = coefficients.get(s34, 0) - 1

        for seg in coefficients:
            if seg not in self.col_id:
                self.add_col(seg)

        row = [0] * len(self.header)
        for seg, coeff in coefficients.items():
            if coeff != 0:
                row[self.col_id[seg]] = coeff

        self.add_row(row, cong)


## TODO: fix the area to the right config
class AreaTable(Table):    
    # each column represents a triangle formed by the two points and the origin(denoted by ZERO_ZERO)
    def __init__(self, header: List[frozenset[Point]]):
        super().__init__(header)

    def add_eqarea(self, eqarea: EqArea):
        p1, p2, p3, p4, p5, p6 = eqarea.points

        coefficients = {}
        seg1 = frozenset([p1, p2])
        seg2 = frozenset([p1, p3])
        seg3 = frozenset([p2, p3])
        seg4 = frozenset([p4, p5])
        seg5 = frozenset([p4, p6])
        seg6 = frozenset([p5, p6])

        if (is_sameclock(p1, p2, p3, p4, p5, p6)):
            coefficients[seg1] = coefficients.get(seg1, 0) + 1
            coefficients[seg2] = coefficients.get(seg2, 0) + 1
            coefficients[seg3] = coefficients.get(seg3, 0) + 1
            coefficients[seg4] = coefficients.get(seg4, 0) - 1
            coefficients[seg5] = coefficients.get(seg5, 0) - 1
            coefficients[seg6] = coefficients.get(seg6, 0) - 1
        else:
            coefficients[seg1] = coefficients.get(seg1, 0) + 1
            coefficients[seg2] = coefficients.get(seg2, 0) + 1
            coefficients[seg3] = coefficients.get(seg3, 0) + 1
            coefficients[seg4] = coefficients.get(seg4, 0) + 1
            coefficients[seg5] = coefficients.get(seg5, 0) + 1
            coefficients[seg6] = coefficients.get(seg6, 0) + 1

        for seg in coefficients:
            if seg not in self.col_id:
                self.add_col(seg)
        
        row = [0] * len(self.header)
        for seg, coeff in coefficients.items():
            if coeff != 0:
                row[self.col_id[seg]] = coeff
        
        self.add_row(row, eqarea)