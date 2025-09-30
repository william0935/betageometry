## Create the Problem data structure
from relations import *

RELATION_TYPES = ["cong", "eqangle", "sameclock", "para", "perp", "col", "cyclic",
                  "simtri1", "simtri2", "contri1", "contri2", "midp"]

class Problem:
    def __init__(self, name: str, points: List[Point],
                 assumptions: Optional[List[RelationNode]] = None,
                 goals: Optional[List[RelationNode]] = None):
        self.name = name
        self.points = points
        self.assumptions = assumptions if assumptions is not None else []
        self.goals = goals if goals is not None else []
        self.solved = False
        self.relations = {r : [] for r in RELATION_TYPES}
        self.index_counter = 1
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
        return f"Problem {self.name}:\n{points_str}\n{assumptions_str}{relations_str}{goals_str}"

    def is_solved(self):
        return self.solved
    
    def add_relation(self, relation: RelationNode):
        if not relation.name in RELATION_TYPES:
            raise ValueError(f"Relation type {relation.name} not recognized.")
        
        for existing in self.relations[relation.name]:
            if existing.relation == relation.relation:
                return
            
        if relation.index is None:
            relation.add_index(self.index_counter)
            self.index_counter += 1

        self.relations[relation.name].append(relation)

        if relation.equivalent:
            for eq in relation.equivalent:
                self.add_relation(eq)

        if relation in self.goals:
            self.goals.remove(relation)
            if not self.goals:
                self.solved = True