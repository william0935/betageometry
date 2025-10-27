import math
from relations import *
import numpy as np
from typing import List, Tuple, Dict, Any, Optional

# Special column key: represents a fixed +90 degrees constant in column 0
CONST_90 = "90_degrees"

class Table:
    def __init__(self, header: List[Any]):
        self.header = list(header)
        self.col_id: Dict[Any, int] = {col: i for i, col in enumerate(self.header)}
        self.rows: Dict[RelationNode, List[np.ndarray]] = {}
        self.relations: List[RelationNode] = []  # Track all relations in order

    def add_col(self, col_name: Any):
        if col_name in self.col_id:
            return
        self.header.append(col_name)
        self.col_id[col_name] = len(self.header) - 1
        # Extend all existing rows with a zero entry
        for relation_node, row_list in self.rows.items():
            for i in range(len(row_list)):
                row_list[i] = np.append(row_list[i], 0)

    # def is_spanned(self, row: List[Any], tol: float = 1e-10) -> bool:
    #     if not self.rows:
    #         return False
    #     r = np.asarray(row, dtype=float)
    #     # Collect all rows from all relations
    #     all_rows = []
    #     for row_list in self.rows.values():
    #         all_rows.extend(row_list)
    #     R = np.vstack(all_rows)  # shape (k, n)
    #     # Check if r is in the row space of R: find w such that R^T w ≈ r^T
    #     A = R.T  # shape (n, k)
    #     w, _, _, _ = np.linalg.lstsq(A, r, rcond=None)
    #     residual = np.linalg.norm(A @ w - r, ord=2)
    #     return residual <= tol
    
    def is_spanned(self, row: List[Any], tol: float = 1e-10) -> Tuple[bool, List[RelationNode]]:
        """Check if a row is in the span of existing rows using Gaussian elimination without division.
        
        Returns:
            Tuple[bool, List[RelationNode]]: (is_spanned, list of relation keys whose rows were used)
        """
        if not self.rows:
            return False, []
        
        # Collect all rows and track which relation each row belongs to
        all_rows = []
        row_to_relation = []  # Maps each row index to its relation key
        for relation_key, row_list in self.rows.items():
            for r in row_list:
                all_rows.append(r)
                row_to_relation.append(relation_key)
        
        r = np.asarray(row, dtype=int)
        R = np.vstack([np.asarray(existing_row, dtype=int) for existing_row in all_rows])
        augmented = np.vstack([R, r])
        n_rows, n_cols = augmented.shape
        
        # Track which rows participated in eliminating the last row
        # A row "participates" if it was used as a pivot or interacted with the last row
        used_row_indices = set()
        
        pivot_row = 0
        for col in range(n_cols):
            found_pivot = False
            pivot_row_idx = None
            for row_idx in range(pivot_row, n_rows):
                if augmented[row_idx, col] != 0:
                    # Swap rows if needed
                    if row_idx != pivot_row:
                        augmented[[pivot_row, row_idx]] = augmented[[row_idx, pivot_row]]
                        # Also swap in the tracking list
                        if pivot_row < len(row_to_relation) and row_idx < len(row_to_relation):
                            row_to_relation[pivot_row], row_to_relation[row_idx] = row_to_relation[row_idx], row_to_relation[pivot_row]
                    found_pivot = True
                    pivot_row_idx = pivot_row
                    break
            
            if not found_pivot:
                continue
            
            # Eliminate all other rows using this pivot (no division)
            pivot_val = augmented[pivot_row, col]
            last_row_idx = n_rows - 1
            
            for row_idx in range(n_rows):
                if row_idx == pivot_row:
                    continue
                if augmented[row_idx, col] == 0:
                    continue
                
                # If we're eliminating the last row (the test row), record which pivot was used
                if row_idx == last_row_idx and pivot_row < len(row_to_relation):
                    used_row_indices.add(pivot_row)
                
                # Eliminate: row = row * pivot_val - augmented[row, col] * pivot_row
                row_multiplier = augmented[row_idx, col]
                augmented[row_idx] = augmented[row_idx] * pivot_val - row_multiplier * augmented[pivot_row]
            
            pivot_row += 1
            if pivot_row >= n_rows - 1:
                break
        
        # Check if the last row (which was r) became all zeros
        last_row = augmented[-1]
        is_spanned_result = np.allclose(last_row, 0, atol=tol)
        
        # Collect the relation keys that were used
        used_relations = []
        if is_spanned_result:
            seen_relation_ids = set()
            for row_idx in sorted(used_row_indices):
                if row_idx < len(row_to_relation):
                    relation = row_to_relation[row_idx]
                    relation_id = id(relation)
                    if relation_id not in seen_relation_ids:
                        seen_relation_ids.add(relation_id)
                        used_relations.append(relation)
        
        return is_spanned_result, used_relations

    def add_row(self, row: List[Any], relation: RelationNode):
        if len(row) != len(self.header):
            raise ValueError("Row length does not match header length.")
        is_spanned, _ = self.is_spanned(row)
        if not is_spanned:
            new_row = np.asarray(row, dtype=int)
            if relation not in self.rows:
                self.rows[relation] = []
            self.rows[relation].append(new_row)
            self.relations.append(relation)

    def table_length(self) -> int:
        return len(self.header)

class AngleTable(Table):
    def __init__(self, header: List[Any]):
        # Ensure the first column is the 90-degree constant
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
        # perp A B C D => AB - CD - 90 = 0
        p1, p2, p3, p4 = perpendicular.points
        seg1 = frozenset({p1, p2})
        seg2 = frozenset({p3, p4})

        angle1 = math.atan2(p2.y - p1.y, p2.x - p1.x)
        angle2 = math.atan2(p4.y - p3.y, p4.x - p3.x)
        angle1 = (angle1 + math.pi) % math.pi
        angle2 = (angle2 + math.pi) % math.pi

        # print(p1, p2, p3, p4, angle1 / math.pi * 180, angle2 / math.pi * 180)

        if angle1 < angle2:
            seg1, seg2 = seg2, seg1 # ensure consistent ordering, not sure if correct

        coefficients: Dict[Any, int] = {}
        coefficients[seg1] = coefficients.get(seg1, 0) + 1
        coefficients[seg2] = coefficients.get(seg2, 0) - 1
        coefficients[CONST_90] = coefficients.get(CONST_90, 0) - 1  # subtract 90 constant

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