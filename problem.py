## Create the Problem data structure
from relations import *
from typing import List, Optional, Tuple
from itertools import combinations, permutations
import numpy as np

RELATION_TYPES = ["cong", "eqangle", "para", "perp", "col", "cyclic", "eqratio",
                  "simtri1", "simtri2", "contri1", "contri2", "midp", "circle", "eqarea"]

class Problem:
    def __init__(self, name: str, points: List[Point],
                 assumptions: Optional[List[RelationNode]] = None,
                 goals: Optional[List[RelationNode]] = None):
        self.name = name
        self.points = points
        self.assumptions = assumptions if assumptions is not None else []
        self.goals = goals if goals is not None else []
        self.remaining_goals = goals.copy() if goals is not None else []
        self.solved = False
        self.relations = {r : [] for r in RELATION_TYPES}
        self.index_counter = 1
        self.deduction_steps = []
        self.similar_triangle_pairs = self.find_similar_triangle_pairs()
        self.cyclic_quads = self.find_cyclic_quads()
        print(self.cyclic_quads)
        for r in self.assumptions:
            r.add_index(self.index_counter)
            self.index_counter += 1

        for r in self.assumptions:
            self.add_relation(r)

        # TODO: implement pairs when needed for AR
        # self.pairs = {(p1, p2) for i, p1 in enumerate(points) for p2 in points[i+1:]}

    def __repr__(self):
        points_str = "Points: " + ", ".join(repr(p) for p in sorted(self.points, key=lambda p: p.name))
        assumptions_str = "Assumptions:\n" + "".join(f"{r}\n" for r in self.assumptions)
        relations_str = "Known Relations:\n"
        relation_list = []
        for r_type in RELATION_TYPES:
            relation_list.extend(self.relations[r_type])
        relations_str += "".join(f"{r}\n" for r in sorted(relation_list, key=lambda r: r.index))
        goals_str = "Goals (solved):\n" if self.solved else "Goals (unsolved):\n"
        goals_str += "".join(f"{g.representation}\n" for g in self.goals)
        procedure_str = ""
        if self.solved:
            procedure_str = self.trace_back()

        # test_output prints all known relations
        test_output = f"Problem {self.name}:\n{points_str}\n{assumptions_str}{relations_str}{goals_str}{procedure_str}"

        # normal_output only prints deduction steps involved in the proof
        normal_output = f"Problem {self.name}:\n{points_str}\n{assumptions_str}{goals_str}{procedure_str}"
        return test_output
    
    def is_solved(self):
        return self.solved

    def add_relation(self, relation: RelationNode) -> Optional[List[RelationNode]]:
        new_relations = []
        if not relation.name in RELATION_TYPES:
            raise ValueError(f"Relation type {relation.name} not recognized.")
        
        if self.is_relation_trivial(relation):
            return None

        for existing in self.relations[relation.name]:
            if existing.relation == relation.relation:
                return None

        if relation.index is None:
            relation.add_index(self.index_counter)
            self.index_counter += 1

        self.relations[relation.name].append(relation)
        new_relations.append(relation)

        if relation.equivalent:
            for eq in relation.equivalent:
                new_rel = self.add_relation(eq)
                if new_rel:
                    new_relations.extend(new_rel)

        for goal in self.remaining_goals:
            if goal.name == relation.name and goal.relation == relation.relation:
                self.remaining_goals.remove(goal)
                self.deduction_steps.append(relation)
                if not self.remaining_goals:
                    self.solved = True
                break

        return new_relations

    def is_relation_trivial(self, relation: RelationNode) -> bool:
        if relation.name == "cong":
            p1, p2, p3, p4 = relation.points
            return set({p1, p2}) == set({p3, p4})
        elif relation.name == "eqangle":
            p1, p2, p3, p4, p5, p6 = relation.points
            return (frozenset({p1, p2}), frozenset({p2, p3})) == (frozenset({p4, p5}), frozenset({p5, p6}))
        elif relation.name == "contri1" or relation.name == "contri2":
            p1, p2, p3, p4, p5, p6 = relation.points
            return (p1, p2, p3) == (p4, p5, p6)
        elif relation.name == "simtri1" or relation.name == "simtri2":
            p1, p2, p3, p4, p5, p6 = relation.points
            return set((p1, p2, p3)) == set((p4, p5, p6)) or \
                   any(relation.relation == contri.relation for contri in self.relations["contri1"] + self.relations["contri2"])
        return False
    
    def trace_back(self) -> str:
        if not self.solved:
            return ""
        
        next_steps = self.deduction_steps.copy()
        while next_steps:
            step = next_steps.pop(0)
            for parent in step.parents:
                if parent not in next_steps:
                    next_steps.append(parent)
            
            if step not in self.deduction_steps:
                self.deduction_steps.append(step)

        steps_str = "Deduction Steps (in order of discovery):\n"
        step_cnt = 1
        for step in sorted(self.deduction_steps, key=lambda r: r.index):
            # reassign indices to output only the steps used in the proof
            step.add_index(step_cnt)
            step_cnt += 1
            steps_str += f"{step}\n"

        return steps_str

    def find_similar_triangle_pairs(self, tol=1e-5) -> List[Tuple[Point, Point, Point, Point, Point, Point]]:
        "find all possible similar triangle pairs among the points according to the diagram"
        points = self.points
        pairs = []
        for p1, p2, p3 in combinations(points, 3):
            for p4, p5, p6 in permutations(points, 3):
                if p1 == p4 and p2 == p5 and p3 == p6:
                    continue
                if abs(self.angle_value(p1, p2, p3) - self.angle_value(p4, p5, p6)) < tol and \
                   abs(self.angle_value(p2, p3, p1) - self.angle_value(p5, p6, p4)) < tol and \
                   abs(self.angle_value(p3, p1, p2) - self.angle_value(p6, p4, p5)) < tol:
                    pairs.append((p1, p2, p3, p4, p5, p6))
        return pairs
    
    def find_cyclic_quads(self, tol=1e-5) -> List[Tuple[Point, Point, Point, Point]]:
        "find all possible cyclic quadrilaterals among the points according to the diagram"
        points = self.points
        quads = []
        for p1, p2 in combinations(points, 2):
            for p3, p4 in permutations(points, 2):
                if len({p1, p2, p3, p4}) < 4:
                    continue
                angle1 = self.angle_value(p1, p3, p2)
                angle2 = self.angle_value(p1, p4, p2)
                if abs(angle1 + angle2 - 180) < tol or abs(angle1 - angle2) < tol:
                    quads.append((p1, p2, p3, p4))
        return quads

    def angle_value(self, a: Point, b: Point, c: Point) -> float:
        "compute the angle value of angle ABC in degrees"
        A = np.array([a.x, a.y])
        B = np.array([b.x, b.y])
        C = np.array([c.x, c.y])
        BA = A - B
        BC = C - B
        cos_angle = np.dot(BA, BC) / (np.linalg.norm(BA) * np.linalg.norm(BC))
        angle_rad = np.arccos(np.clip(cos_angle, -1.0, 1.0))
        angle_deg = np.degrees(angle_rad)
        # normalize angle degree to [0, 180]
        if angle_deg > 180:
            angle_deg = 360 - angle_deg
        return angle_deg