from relations import *
import numpy as np
from typing import List, Tuple, Dict, Any, Optional

class Table:
    def __init__(self, header: List[Any]):
        self.header = header
        self.col_id: Dict[Any, int] = {col: i for i, col in enumerate(header)}
        self.relations: Dict[int, RelationNode] = {}
        self.rows: List[np.ndarray] = []  # List of row vectors in RREF form
        self.pivot_cols: List[int] = []  # Track pivot columns for each row

    def add_row(self, row: List[Any], relation: RelationNode):
        if len(row) != len(self.header):
            raise ValueError("Row length does not match header length.")
        
        new_row = np.array(row, dtype=float)
        self._reduce_row_inplace(new_row)
        pivot_col = self._find_pivot(new_row)
        
        if pivot_col is not None:
            new_row /= new_row[pivot_col]
            for i, existing_row in enumerate(self.rows):
                if abs(existing_row[pivot_col]) > 1e-10:
                    existing_row -= existing_row[pivot_col] * new_row
            
            insert_pos = self._find_insert_position(pivot_col)
            self.rows.insert(insert_pos, new_row)
            self.pivot_cols.insert(insert_pos, pivot_col)
            
            new_relations = {}
            for idx, rel in self.relations.items():
                if idx >= insert_pos:
                    new_relations[idx + 1] = rel
                else:
                    new_relations[idx] = rel
            new_relations[insert_pos] = relation
            self.relations = new_relations

    def add_col(self, col_name: Any):
        self.header.append(col_name)

        for row in self.rows:
            row.resize(len(self.header), refcheck=False)
            row[-1] = 0
        self.col_id[col_name] = len(self.header) - 1

    def _reduce_row_inplace(self, row: np.ndarray):
        for i, pivot_col in enumerate(self.pivot_cols):
            if abs(row[pivot_col]) > 1e-10:
                row -= row[pivot_col] * self.rows[i]

    def _find_pivot(self, row: np.ndarray) -> Optional[int]:
        for i, val in enumerate(row):
            if abs(val) > 1e-10:
                return i
        return None

    def _find_insert_position(self, pivot_col: int) -> int:
        left, right = 0, len(self.pivot_cols)
        while left < right:
            mid = (left + right) // 2
            if self.pivot_cols[mid] < pivot_col:
                left = mid + 1
            else:
                right = mid
        return left

    def is_spanned(self, row: List[Any]) -> bool:
        if not self.rows:
            return False
        
        test_row = np.array(row, dtype=float)
        self._reduce_row_inplace(test_row)
        return np.allclose(test_row, 0, atol=1e-10)
    
    def table_length(self) -> int:
        return len(self.header)

class AngleTable(Table):
    def __init__(self, header: List[frozenset[Point, Point]]):
        super().__init__(header)
    
    def add_eqangle(self, eqangle: EqualAngle):
        p1, p2, p3, p4, p5, p6 = eqangle.points
        seg1_1 = frozenset({p1, p2})
        seg1_2 = frozenset({p2, p3})
        seg2_1 = frozenset({p4, p5})
        seg2_2 = frozenset({p5, p6})
        
        # Using coefficients handles the case of duplicate points
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

        if not self.is_spanned(row):
            self.add_row(row, eqangle)
    
    def add_collinear(self, collinear : Collinear):
        p1, p2, p3 = collinear.points
        seg1 = frozenset({p1, p2})
        seg2 = frozenset({p2, p3})

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

        if not self.is_spanned(row):
            self.add_row(row, collinear)

class RatioTable(Table):
    def __init__(self, header: List[frozenset[Point, Point]]):
        super().__init__(header)
    
    def add_eqratio(self, eqratio: EqualRatio):
        p1, p2, p3, p4, p5, p6, p7, p8 = eqratio.points
        row = [0] * len(self.header)
        if frozenset({p1, p2}) not in self.col_id:
            self.add_col(frozenset({p1, p2}))
        if frozenset({p3, p4}) not in self.col_id:
            self.add_col(frozenset({p3, p4}))
        if frozenset({p5, p6}) not in self.col_id:
            self.add_col(frozenset({p5, p6}))
        if frozenset({p7, p8}) not in self.col_id:
            self.add_col(frozenset({p7, p8}))

        row[self.col_id[frozenset({p1, p2})]] = 1
        row[self.col_id[frozenset({p3, p4})]] = -1
        row[self.col_id[frozenset({p5, p6})]] = -1
        row[self.col_id[frozenset({p7, p8})]] = 1

        if not self.is_spanned(row):
            self.add_row(row, eqratio)

    def add_cong(self, cong: Congruent):
        p1, p2, p3, p4 = cong.points
        row = [0] * len(self.header)
        if frozenset({p1, p2}) not in self.col_id:
            self.add_col(frozenset({p1, p2}))
        if frozenset({p3, p4}) not in self.col_id:
            self.add_col(frozenset({p3, p4}))

        row[self.col_id[frozenset({p1, p2})]] = 1
        row[self.col_id[frozenset({p3, p4})]] = -1

        if not self.is_spanned(row):
            self.add_row(row, cong)
