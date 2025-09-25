## Create the Problem data structure
from relations import *

RELATION_TYPES = ["cong", "eqangle", "sameclock", "para", "perp", "col", "cyclic",
                  "simtri1", "simtri2", "contri1", "contri2", "midp"]

class Problem:
    def __init__(self, name: str, points: set[Point],
                 assumptions: Optional[List[RelationNode]] = None,
                 goals: Optional[List[RelationNode]] = None):
        self.name = name
        self.points = points
        self.assumptions = assumptions if assumptions is not None else []
        self.goals = goals if goals is not None else []
        self.solved = False
        self.relations = {r : [] for r in RELATION_TYPES}
        for r in self.given_relations:
            self.add_relation(r)

    def __repr__(self):
        points_str = "Points: " + ", ".join(repr(p) for p in sorted(self.points, key=lambda p: p.name))
        given_relations_str = "Assumptions:\n" + "".join(f"[{i+1}]: {r}\n" for i, r in enumerate(self.given_relations))
        relations_id = 1
        relations_str = "Known Relations:\n"
        for r_type in RELATION_TYPES:
            relations_str += "".join(f"[{relations_id+i}]: {r}\n" for i, r in enumerate(self.relations[r_type]))
            relations_id += len(self.relations[r_type])
        goals_str = "Goals (solved):\n" if self.solved else "Goals (unsolved):\n"
        goals_str += "".join(f"[{i+1}]: {g}\n" for i, g in enumerate(self.goals))
        return f"Problem {self.name}:\n{points_str}\n{given_relations_str}{relations_str}{goals_str}"

    def is_solved(self):
        return self.solved
    
    def add_relation(self, relation: RelationNode):
        if not relation.name in RELATION_TYPES:
            raise ValueError(f"Relation type {relation.name} not recognized.")
        
        for existing in self.relations[relation.name]:
            if existing.relation == relation.relation:
                return
            
        self.relations[relation.name].append(relation)
        for eq in relation.equivalent:
            self.add_relation(eq)
        
        if relation in self.goals:
            self.goals.remove(relation)
            if not self.goals:
                self.solved = True