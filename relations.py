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
    def __init__(self, name: str, index: int, relation, parents: Optional[set["RelationNode"]] = None, rule: Optional[str] = None,
                 representation: str = "", equivalent: Optional[List["RelationNode"]] = None):
        self.name = name
        self.index = index
        self.relation = relation
        self.parents = parents if parents is not None else set()
        self.rule = rule if rule is not None else ""
        self.representation = representation
        self.equivalent = equivalent if equivalent is not None else []
        self.points = []

    def add_index(self, index: int):
        self.index = index
        return self
    
    # TODO: improve representation
    def __repr__(self):
        # return f"[{self.index}] {self.representation} from [{' ,'.join(str(p.index) for p in self.parents)}] via {self.rule}" if self.rule else f"[{self.index}] {self.representation}"
        repr_str = f"[{self.index}] {self.representation}"
        if self.parents:
            repr_str += f" from [{', '.join(str(p.index) for p in self.parents)}]"
        if self.rule:
            repr_str += f" via {self.rule}"
        return repr_str


class Congruent(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point, p4: Point,
                 index: Optional[int] = None,
                 parents: Optional[set[RelationNode]] = None,
                 rule: Optional[str] = None):
        super().__init__(
            name="cong",
            index=index,
            relation=frozenset({frozenset({p1, p2}), frozenset({p3, p4})}),
            representation=f"cong {p1} {p2} {p3} {p4}",
            parents=parents,
            rule=rule
        )

        self.points = [p1, p2, p3, p4]


class EqualAngle(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point,
                 p4: Point, p5: Point, p6: Point,
                 index: Optional[int] = None,
                 parents: Optional[set[RelationNode]] = None,
                 rule: Optional[str] = None):
        super().__init__(
            name="eqangle",
            index=index,
            relation=frozenset({(frozenset({p1, p2}), frozenset({p2, p3})), 
                                (frozenset({p4, p5}), frozenset({p5, p6}))}),
            representation=f"eqangle {p1} {p2} {p3} {p4} {p5} {p6}",
            parents=parents,
            rule=rule
        )

        self.points = [p1, p2, p3, p4, p5, p6]


# class Sameclock(RelationNode):
#     def __init__(self, p1: Point, p2: Point, p3: Point,
#                  p4: Point, p5: Point, p6: Point,
#                  index: Optional[int] = None,
#                  parents: Optional[set[RelationNode]] = None,
#                  rule: Optional[str] = None):
#         super().__init__(
#             name="sameclock",
#             index=index,
#             relation=frozenset({(p1, p2, p3), (p4, p5, p6)}),
#             representation=f"sameclock {p1} {p2} {p3} {p4} {p5} {p6}",
#             parents=parents,
#             rule=rule
#         )


class Parallel(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point, p4: Point,
                 index: Optional[int] = None,
                 parents: Optional[set[RelationNode]] = None,
                 rule: Optional[str] = None):
        super().__init__(
            name="para",
            index=index,
            relation=frozenset({frozenset({p1, p2}), frozenset({p3, p4})}),
            representation=f"para {p1} {p2} {p3} {p4}",
            parents=parents,
            rule=rule,
            equivalent=[
                EqualAngle(p1, p2, p3, p4, p3, p2, parents=[self]),
                EqualAngle(p2, p1, p3, p4, p3, p1, parents=[self]),
                EqualAngle(p2, p1, p4, p3, p4, p1, parents=[self]),
                EqualAngle(p1, p2, p4, p3, p4, p2, parents=[self])
            ]
        )

        self.points = [p1, p2, p3, p4]


class Perpendicular(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point, p4: Point,
                 index: Optional[int] = None,
                 parents: Optional[set[RelationNode]] = None,
                 rule: Optional[str] = None):
        super().__init__(
            name="perp",
            index=index,
            relation=frozenset({frozenset({p1, p2}), frozenset({p3, p4})}),
            representation=f"perp {p1} {p2} {p3} {p4}",
            parents=parents,
            rule=rule
        )

        self.points = [p1, p2, p3, p4]


class Collinear(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point,
                 index: Optional[int] = None,
                 parents: Optional[set[RelationNode]] = None,
                 rule: Optional[str] = None):
        super().__init__(
            name="col",
            index=index,
            relation=frozenset({p1, p2, p3}),
            representation=f"col {p1} {p2} {p3}",
            parents=parents,
            rule=rule
        )

        self.points = [p1, p2, p3]


class Cyclic(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point, p4: Point,
                 index: Optional[int] = None,
                 parents: Optional[set[RelationNode]] = None,
                 rule: Optional[str] = None):
        super().__init__(
            name="cyclic",
            index=index,
            relation=frozenset({p1, p2, p3, p4}),
            representation=f"cyclic {p1} {p2} {p3} {p4}",
            parents=parents,
            rule=rule,
            equivalent=[
                EqualAngle(p1, p2, p3, p1, p4, p3, parents=[self]),
                EqualAngle(p2, p3, p4, p2, p1, p4, parents=[self]),
                EqualAngle(p1, p2, p4, p1, p3, p4, parents=[self]),
                EqualAngle(p2, p1, p3, p2, p4, p3, parents=[self])
            ]
        )

        self.points = [p1, p2, p3, p4]


class SimilarTriangle1(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point,
                 p4: Point, p5: Point, p6: Point,
                 index: Optional[int] = None,
                 parents: Optional[set[RelationNode]] = None,
                 rule: Optional[str] = None):
        super().__init__(
            name="simtri1",
            index=index,
            relation=frozenset({
                frozenset({p1, p4}),
                frozenset({p2, p5}),
                frozenset({p3, p6})
            }),
            representation=f"simtri1 {p1} {p2} {p3} {p4} {p5} {p6}",
            parents=parents,
            rule=rule,
            equivalent=[
                EqualAngle(p1, p2, p3, p4, p5, p6, parents=[self]),
                EqualAngle(p2, p3, p1, p5, p6, p4, parents=[self]),
                EqualAngle(p3, p1, p2, p6, p4, p5, parents=[self]),
                EqualRatio(p1, p2, p2, p3, p4, p5, p5, p6, parents=[self]),
                EqualRatio(p2, p3, p3, p1, p5, p6, p6, p4, parents=[self]),
                EqualRatio(p3, p1, p1, p2, p6, p4, p4, p5, parents=[self]),
                # Sameclock(p1, p2, p3, p4, p5, p6, parents=[self])
            ]
        )

        self.points = [p1, p2, p3, p4, p5, p6]


class SimilarTriangle2(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point,
                 p4: Point, p5: Point, p6: Point,
                 index: Optional[int] = None,
                 parents: Optional[set[RelationNode]] = None,
                 rule: Optional[str] = None):
        super().__init__(
            name="simtri2",
            index=index,
            relation=frozenset({
                frozenset({p1, p4}),
                frozenset({p2, p5}),
                frozenset({p3, p6})
            }),
            representation=f"simtri2 {p1} {p2} {p3} {p4} {p5} {p6}",
            parents=parents,
            rule=rule,
            equivalent=[
                EqualAngle(p1, p2, p3, p6, p5, p4, parents=[self]),
                EqualAngle(p2, p3, p1, p4, p6, p5, parents=[self]),
                EqualAngle(p3, p1, p2, p5, p4, p6, parents=[self]),
                EqualRatio(p1, p2, p2, p3, p4, p5, p5, p6, parents=[self]),
                EqualRatio(p2, p3, p3, p1, p5, p6, p6, p4, parents=[self]),
                EqualRatio(p3, p1, p1, p2, p6, p4, p4, p5, parents=[self]),
                # Sameclock(p1, p2, p3, p6, p5, p4, parents=[self])
            ]
        )

        self.points = [p1, p2, p3, p4, p5, p6]


class CongruentTriangle1(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point,
                 p4: Point, p5: Point, p6: Point,
                 index: Optional[int] = None,
                 parents: Optional[set[RelationNode]] = None,
                 rule: Optional[str] = None):
        super().__init__(
            name="contri1",
            index=index,
            relation=frozenset({
                frozenset({p1, p4}),
                frozenset({p2, p5}),
                frozenset({p3, p6})
            }),
            representation=f"contri1 {p1} {p2} {p3} {p4} {p5} {p6}",
            parents=parents,
            rule=rule,
            equivalent=[
                Congruent(p1, p2, p4, p5, parents=[self]),
                Congruent(p2, p3, p5, p6, parents=[self]),
                Congruent(p3, p1, p6, p4, parents=[self]),
                EqualAngle(p1, p2, p3, p4, p5, p6, parents=[self]),
                EqualAngle(p2, p3, p1, p5, p6, p4, parents=[self]),
                EqualAngle(p3, p1, p2, p6, p4, p5, parents=[self]),
                # Sameclock(p1, p2, p3, p4, p5, p6, parents=[self])
            ]
        )

        self.points = [p1, p2, p3, p4, p5, p6]


class CongruentTriangle2(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point,
                 p4: Point, p5: Point, p6: Point,
                 index: Optional[int] = None,
                 parents: Optional[set[RelationNode]] = None,
                 rule: Optional[str] = None):
        super().__init__(
            name="contri2",
            index=index,
            relation=frozenset({
                frozenset({p1, p4}),
                frozenset({p2, p5}),
                frozenset({p3, p6})
            }),
            representation=f"contri2 {p1} {p2} {p3} {p4} {p5} {p6}",
            parents=parents,
            rule=rule,
            equivalent=[
                Congruent(p1, p2, p4, p5, parents=[self]),
                Congruent(p2, p3, p5, p6, parents=[self]),
                Congruent(p3, p1, p6, p4, parents=[self]),
                EqualAngle(p1, p2, p3, p6, p5, p4, parents=[self]),
                EqualAngle(p2, p3, p1, p4, p6, p5, parents=[self]),
                EqualAngle(p3, p1, p2, p5, p4, p6, parents=[self]),
                # Sameclock(p1, p2, p3, p6, p5, p4, parents=[self])
            ]
        )

        self.points = [p1, p2, p3, p4, p5, p6]


class Midpoint(RelationNode):
    def __init__(self, mid: Point, p1: Point, p2: Point,
                 index: Optional[int] = None,
                 parents: Optional[set[RelationNode]] = None,
                 rule: Optional[str] = None):
        super().__init__(
            name="midp",
            index=index,
            relation=(mid, frozenset({p1, p2})),
            representation=f"midp {mid} {p1} {p2}",
            parents=parents,
            rule=rule,
            equivalent=[
                Congruent(mid, p1, mid, p2, parents=[self]),
                Collinear(mid, p1, p2, parents=[self])
            ]
        )

        self.points = [mid, p1, p2]


class EqualRatio(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point, p4: Point,
                 p5: Point, p6: Point, p7: Point, p8: Point,
                 index: Optional[int] = None,
                 parents: Optional[set[RelationNode]] = None,
                 rule: Optional[str] = None):
        super().__init__(
            name="eqratio",
            index=index,
            relation=((frozenset({p1, p2}), frozenset({p3, p4})),
                      (frozenset({p5, p6}), frozenset({p7, p8}))),
            representation=f"eqratio {p1} {p2} {p3} {p4} {p5} {p6} {p7} {p8}",
            parents=parents,
            rule=rule
        )

        self.points = [p1, p2, p3, p4, p5, p6, p7, p8]


class Circle(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point, p4: Point,
                 index: Optional[int] = None,
                 parents: Optional[set[RelationNode]] = None,
                 rule: Optional[str] = None):
        super().__init__(
            name="circle",
            index=index,
            relation=(p1, frozenset({p2, p3, p4})),
            representation=f"circle {p1} {p2} {p3} {p4}",
            parents=parents,
            rule=rule,
            equivalent=[
                Congruent(p1, p2, p1, p3, parents=[self]),
                Congruent(p1, p2, p1, p4, parents=[self]),
                Congruent(p1, p3, p1, p4, parents=[self])
            ]
        )

        self.points = [p1, p2, p3, p4]


if __name__ == "__main__":
    pass
