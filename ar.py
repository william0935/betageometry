from relations import *
import numpy as np
from typing import List, Tuple, Dict, Any

class Table:
    def __init__(self, header: List[Any]):
        self.header = header
        self.col_id: Dict[Any, int] = {col: i for i, col in enumerate(header)}
        self.rows = np.array([]).reshape(0, len(header))
        self.relations: Dict[int, RelationNode] = {}

    def add_row(self, row: List[Any], relation: RelationNode):
        if len(row) != len(self.header):
            raise ValueError("Row length does not match header length.")
        self.rows = np.vstack([self.rows, row])
        self.relations[self.rows.shape[0] - 1] = relation

    def add_col(self, col_name: Any):
        self.header.append(col_name)
        self.rows = np.hstack([self.rows, np.zeros((self.rows.shape[0], 1))])
        self.col_id.update({col_name: len(self.header) - 1})

    def _rref(self, matrix: np.ndarray) -> np.ndarray:
        A = matrix.copy().astype(float)
        rows, cols = A.shape
        r = 0
        for c in range(cols):
            if r >= rows:
                break
            pivot = np.argmax(np.abs(A[r:rows, c])) + r
            if A[pivot, c] == 0:
                continue
            A[[r, pivot]] = A[[pivot, r]]
            A[r] = A[r] / A[r, c]
            for i in range(rows):
                if i != r:
                    A[i] = A[i] - A[i, c] * A[r]
            r += 1
        return A
    
    def is_spanned(self, row: List[Any]) -> bool:
        if self.rows.shape[0] == 0:
            return False
        
        augmented = np.vstack([self.rows, row])
        rref = self._rref(augmented)
        return np.all(rref[-1, :-1] == 0) and rref[-1, -1] == 0


class AngleTable(Table):
    def __init__(self, header: List[frozenset[Point, Point]]):
        super().__init__(header)
    
    def add_eqangle(self, eqangle: EqualAngle):
        p1, p2, p3, p4, p5, p6 = eqangle.points
        row = [0] * len(self.header)
        if frozenset({p1, p2}) not in self.col_id:
            self.add_col(frozenset({p1, p2}))
        if frozenset({p2, p3}) not in self.col_id:
            self.add_col(frozenset({p2, p3}))
        if frozenset({p4, p5}) not in self.col_id:
            self.add_col(frozenset({p4, p5}))
        if frozenset({p5, p6}) not in self.col_id:
            self.add_col(frozenset({p5, p6}))

        row[self.col_id[frozenset({p1, p2})]] = 1
        row[self.col_id[frozenset({p2, p3})]] = -1
        row[self.col_id[frozenset({p4, p5})]] = -1
        row[self.col_id[frozenset({p5, p6})]] = 1

        if not self.is_spanned(row):
            self.add_row(row, eqangle)

class RatioTable(Table):
    def __init__(self, header: List[frozenset[Point, Point]]):
        super().__init__(header)
    
    def add_eqratio(self, eqratio: EqualRatio):
        p1, p2, p3, p4, p5, p6, p7, p8 = eqratio.points
        row = [0] * len(self.header)
        if (p1, p2) not in self.col_id:
            self.add_col((p1, p2))
        if (p3, p4) not in self.col_id:
            self.add_col((p3, p4))
        if (p5, p6) not in self.col_id:
            self.add_col((p5, p6))
        if (p7, p8) not in self.col_id:
            self.add_col((p7, p8))

        row[self.col_id[(p1, p2)]] = 1
        row[self.col_id[(p3, p4)]] = -1
        row[self.col_id[(p5, p6)]] = -1
        row[self.col_id[(p7, p8)]] = 1

        if not self.is_spanned(row):
            self.add_row(row, eqratio)

    def add_cong(self, cong: Congruent):
        p1, p2, p3, p4 = cong.points
        row = [0] * len(self.header)
        if (p1, p2) not in self.col_id:
            self.add_col((p1, p2))
        if (p3, p4) not in self.col_id:
            self.add_col((p3, p4))

        row[self.col_id[(p1, p2)]] = 1
        row[self.col_id[(p3, p4)]] = -1

        if not self.is_spanned(row):
            self.add_row(row, cong)
