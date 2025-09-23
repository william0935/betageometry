from typing import List, Optional


# Basic structures
class Point:
    def __init__(self, name: str, x: float = None, y: float = None):
        self.name = name
        self.x = x
        self.y = y

    def __repr__(self):
        return self.name


# Reference objects
ZERO_ZERO = Point("_0", 0.0, 0.0)
UNIT_X = Point("_x", 1.0, 0.0)
UNIT_Y = Point("_y", 0.0, 1.0)
MINUS_UNIT_X = Point("_-x", -1.0, 0.0)

ZERO_ANGLE = (frozenset({ZERO_ZERO, UNIT_X}), frozenset({ZERO_ZERO, UNIT_X}))
PERPENDICULAR_ANGLE = (frozenset({ZERO_ZERO, UNIT_X}), frozenset({ZERO_ZERO, UNIT_Y}))


# Base class for relations
class RelationNode:
    def __init__(self, name: str, relation, parents: Optional[List["RelationNode"]] = None,
                 representation: str = "", equivalent: Optional[List["RelationNode"]] = None):
        self.name = name
        self.relation = relation
        self.parents = parents if parents is not None else []
        self.representation = representation
        self.equivalent = equivalent if equivalent is not None else []

    def __repr__(self):
        parent_str = ""
        if self.parents:
            parent_str = "\nfrom:"
            for p in self.parents:
                parent_str += f"\n  {repr(p)}"
        return f"{self.representation}" + parent_str


class Congruent(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point, p4: Point,
                 parents: Optional[List[RelationNode]] = None):
        super().__init__(
            name="cong",
            relation=frozenset({frozenset({p1, p2}), frozenset({p3, p4})}),
            representation=f"cong {p1} {p2} {p3} {p4}",
            parents=parents
        )


class EqualAngle(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point,
                 p4: Point, p5: Point, p6: Point,
                 parents: Optional[List[RelationNode]] = None):
        super().__init__(
            "eqangle",
            relation=frozenset({
                (frozenset({p1, p2}), frozenset({p2, p3})),
                (frozenset({p4, p5}), frozenset({p5, p6}))
            }),
            representation=f"eqangle {p1} {p2} {p3} {p4} {p5} {p6}",
            parents=parents
        )


class Sameclock(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point,
                 p4: Point, p5: Point, p6: Point,
                 parents: Optional[List[RelationNode]] = None):
        super().__init__(
            "sameclock",
            relation=frozenset({(p1, p2, p3), (p4, p5, p6)}),
            representation=f"sameclock {p1} {p2} {p3} {p4} {p5} {p6}",
            parents=parents
        )


class Parallel(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point, p4: Point,
                 parents: Optional[List[RelationNode]] = None):
        super().__init__(
            "para",
            relation=frozenset({(p1, p2), (p3, p4)}),
            representation=f"para {p1} {p2} {p3} {p4}",
            parents=parents
        )
        self.equivalent = [EqualAngle(p1, p2, p3, p4, p3, p2, [self])]


class Perpendicular(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point, p4: Point,
                 parents: Optional[List[RelationNode]] = None):
        super().__init__(
            "perp",
            relation=frozenset({(p1, p2), (p3, p4)}),
            representation=f"perp {p1} {p2} {p3} {p4}",
            parents=parents
        )
        # TODO: define equivalent relations if needed


class Colinear(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point,
                 parents: Optional[List[RelationNode]] = None):
        super().__init__(
            "col",
            relation=frozenset({p1, p2, p3}),
            representation=f"col {p1} {p2} {p3}",
            parents=parents
        )
        self.equivalent = [EqualAngle(p1, p2, p3, p3, p2, p1, [self])]


class Cyclic(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point, p4: Point,
                 parents: Optional[List[RelationNode]] = None):
        super().__init__(
            "cyclic",
            relation=frozenset({p1, p2, p3, p4}),
            representation=f"cyclic {p1} {p2} {p3} {p4}",
            parents=parents
        )
        self.equivalent = [EqualAngle(p1, p2, p3, p1, p4, p3, [self])]


class SimilarTriangle1(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point,
                 p4: Point, p5: Point, p6: Point,
                 parents: Optional[List[RelationNode]] = None):
        super().__init__(
            "simtri1",
            relation=frozenset({
                frozenset({p1, p4}),
                frozenset({p2, p5}),
                frozenset({p3, p6})
            }),
            representation=f"simtri1 {p1} {p2} {p3} {p4} {p5} {p6}",
            parents=parents
        )
        self.equivalent = [
            EqualAngle(p1, p2, p3, p4, p5, p6, [self]),
            EqualAngle(p2, p3, p1, p5, p6, p4, [self]),
            Sameclock(p1, p2, p3, p4, p5, p6, [self])
        ]


class SimilarTriangle2(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point,
                 p4: Point, p5: Point, p6: Point,
                 parents: Optional[List[RelationNode]] = None):
        super().__init__(
            "simtri2",
            relation=frozenset({
                frozenset({p1, p4}),
                frozenset({p2, p5}),
                frozenset({p3, p6})
            }),
            representation=f"simtri2 {p1} {p2} {p3} {p4} {p5} {p6}",
            parents=parents
        )
        self.equivalent = [
            EqualAngle(p1, p2, p3, p6, p5, p4, [self]),
            EqualAngle(p2, p3, p1, p4, p6, p5, [self]),
            Sameclock(p1, p2, p3, p6, p5, p4, [self])
        ]


class CongruentTriangle1(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point,
                 p4: Point, p5: Point, p6: Point,
                 parents: Optional[List[RelationNode]] = None):
        super().__init__(
            "contri1",
            relation=frozenset({
                frozenset({p1, p4}),
                frozenset({p2, p5}),
                frozenset({p3, p6})
            }),
            representation=f"contri1 {p1} {p2} {p3} {p4} {p5} {p6}",
            parents=parents
        )
        self.equivalent = [
            Congruent(p1, p2, p4, p5, [self]),
            Congruent(p2, p3, p5, p6, [self]),
            EqualAngle(p1, p2, p3, p4, p5, p6, [self]),
            Sameclock(p1, p2, p3, p4, p5, p6, [self])
        ]


class CongruentTriangle2(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point,
                 p4: Point, p5: Point, p6: Point,
                 parents: Optional[List[RelationNode]] = None):
        super().__init__(
            "contri2",
            relation=frozenset({
                frozenset({p1, p4}),
                frozenset({p2, p5}),
                frozenset({p3, p6})
            }),
            representation=f"contri2 {p1} {p2} {p3} {p4} {p5} {p6}",
            parents=parents
        )
        self.equivalent = [
            Congruent(p1, p2, p4, p5, [self]),
            Congruent(p2, p3, p5, p6, [self]),
            EqualAngle(p1, p2, p3, p6, p5, p4, [self]),
            Sameclock(p1, p2, p3, p6, p5, p4, [self])
        ]


class Midpoint(RelationNode):
    def __init__(self, mid: Point, p1: Point, p2: Point,
                 parents: Optional[List[RelationNode]] = None):
        super().__init__(
            "midp",
            relation=frozenset({mid, p1, p2}),
            representation=f"midp {mid} {p1} {p2}",
            parents=parents
        )
        self.equivalent = [
            Congruent(mid, p1, mid, p2, [self]),
            Colinear(mid, p1, p2, [self])
        ]


# Problem structure
class Problem:
    def __init__(self, name: str, points: set[Point],
                 relations: Optional[List[RelationNode]] = None,
                 goals: Optional[List[RelationNode]] = None):
        self.name = name
        self.points = points
        self.relations = relations if relations is not None else []
        self.goals = goals if goals is not None else []
        self.solved = False
        self.expand_relations()

    def __repr__(self):
        points_str = ", ".join(repr(p) for p in sorted(self.points, key=lambda p: p.name))
        relations_str = "".join(f"[{i+1}]: {r}\n" for i, r in enumerate(self.relations))
        goals_str = " (solved):\n" if self.solved else "(unsolved):\n"
        goals_str += ":".join(f"[{i+1}]: {g}\n" for i, g in enumerate(self.goals))
        return f"Problem {self.name}:\nPoints: {points_str}\nRelations:\n{relations_str}Goals{goals_str}"

    def expand_relations(self):
        new_relations = []
        for r in self.relations:
            new_relations.extend(r.equivalent)
        self.relations.extend(new_relations)

    def is_solved(self):
        return self.solved


if __name__ == "__main__":
    pass