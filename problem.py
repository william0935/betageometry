## Create the Problem data structure
from relations import *

RELATION_TYPES = ["cong", "eqangle", "para", "perp", "col", "cyclic", "eqratio",
                  "simtri1", "simtri2", "contri1", "contri2", "midp", "circle"]

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
            if goal.relation == relation.relation:
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
            return (p1, p2, p3) == (p4, p5, p6) or \
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