from relations import *

RELATION_TYPES = ["cong", "eqangle", "sameclock", "para", "perp", "col", "cyclic",
                  "simtri1", "simtri2", "contri1", "contri2", "midp"]

class Problem:
    def __init__(self, name: str, points: set[Point],
                 given_relations: Optional[List[RelationNode]] = None,
                 goals: Optional[List[RelationNode]] = None):
        self.name = name
        self.points = points
        self.given_relations = given_relations if given_relations is not None else []
        self.goals = goals if goals is not None else []
        self.solved = False
        self.relations = {r : [] for r in RELATION_TYPES}
        self.get_relations()

    def __repr__(self):
        points_str = "Points: " + ", ".join(repr(p) for p in sorted(self.points, key=lambda p: p.name))
        given_relations_str = "Given Relations:\n" + "".join(f"[{i+1}]: {r}\n" for i, r in enumerate(self.given_relations))
        relations_id = 1
        relations_str = "Relations:\n"
        for r_type in RELATION_TYPES:
            relations_str += "".join(f"[{relations_id+i}]: {r}\n" for i, r in enumerate(self.relations[r_type]))
            relations_id += len(self.relations[r_type])
        goals_str = "Goals (solved):\n" if self.solved else "Goals (unsolved):\n"
        goals_str += "".join(f"[{i+1}]: {g}\n" for i, g in enumerate(self.goals))
        return f"Problem {self.name}:\n{points_str}\n{given_relations_str}{relations_str}{goals_str}"

    def get_relations(self):
        for r in self.given_relations:
            self.relations[r.name].append(r)
            for eq in r.equivalent:
                self.relations[eq.name].append(eq)

    def is_solved(self):
        return self.solved

    
def construct_problem():
    # Define Points
    a = Point("A", 1.0, 2.0)
    b = Point("B", -1.0, -1.0)
    c = Point("C", 2.0, -2.0)
    d = Point("D", 0.5, -1.5)
    e = Point("E", 1.5, 0)
    f = Point("F", 0, 0.5)
    g = Point("G", (a.x + b.x + c.x)/3, (a.y + b.y + c.y)/3)

    # Define Relations
    midp_fab = Midpoint(f, a, b)
    midp_ebc = Midpoint(e, b, c)
    col_agd = Collinear(a, g, d)
    col_bdc = Collinear(b, d, c)
    col_bge = Collinear(b, g, e)
    col_cgf = Collinear(c, g, f)

    goal = Midpoint(g, e, f)

    return Problem(
        name="definition of the centroid",
        points={a, b, c, d, e, f, g},
        given_relations=[midp_fab, midp_ebc, col_agd, col_bdc, col_bge, col_cgf],
        goals=[goal]
    )


def main():
    test_problem = construct_problem()
    print(test_problem)


if __name__ == "__main__":
    main()