import math
from typing import List, Optional


# Basic structures
class Point:
    def __init__(self, name: str, x: float = None, y: float = None):
        self.name = name
        self.x = x
        self.y = y

    def __repr__(self):
        return self.name


# Diagram-level numeric predicates.
#
# These are scale-invariant: they compare a cross/dot product against the product of the
# two lengths, which is sin/cos of the angle between the lines. Comparing gradients
# instead makes the tolerance meaningless for steep lines -- (0,0), (0.001,1),
# (0.002,2.0000001) are collinear but their gradients differ by 1e-4 -- and needs a
# special case for vertical lines that an exact `== 0` test almost never catches on
# coordinates that came out of a file or a construction.

def direction_mod_pi(p1: Point, p2: Point) -> float:
    """Direction of line p1p2 as an angle in [0, pi).

    Normalised so that a direction of pi reports as 0: they are the same line, and
    letting one land at each end of the range makes orientation-sensitive callers
    (see `AngleTable.add_perpendicular`) disagree about the same line.
    """
    a = math.atan2(p2.y - p1.y, p2.x - p1.x) % math.pi
    if a > math.pi - 1e-12:
        a = 0.0
    return a


def angles_close_mod_pi(a: float, b: float, tol: float = 1e-6) -> bool:
    """True when two angles agree modulo pi, including across the 0/pi wrap."""
    d = (a - b) % math.pi
    return min(d, math.pi - d) <= tol


def points_look_collinear(p1: Point, p2: Point, p3: Point, tol: float = 1e-6) -> bool:
    ax, ay = p2.x - p1.x, p2.y - p1.y
    bx, by = p3.x - p1.x, p3.y - p1.y
    na = math.hypot(ax, ay)
    nb = math.hypot(bx, by)
    if na == 0.0 or nb == 0.0:
        return True  # a repeated point is trivially collinear
    return abs(ax * by - ay * bx) <= tol * na * nb


def lines_look_parallel(p1: Point, p2: Point, p3: Point, p4: Point, tol: float = 1e-6) -> bool:
    ax, ay = p2.x - p1.x, p2.y - p1.y
    bx, by = p4.x - p3.x, p4.y - p3.y
    na = math.hypot(ax, ay)
    nb = math.hypot(bx, by)
    if na == 0.0 or nb == 0.0:
        return False
    return abs(ax * by - ay * bx) <= tol * na * nb


def lines_look_perpendicular(p1: Point, p2: Point, p3: Point, p4: Point, tol: float = 1e-6) -> bool:
    ax, ay = p2.x - p1.x, p2.y - p1.y
    bx, by = p4.x - p3.x, p4.y - p3.y
    na = math.hypot(ax, ay)
    nb = math.hypot(bx, by)
    if na == 0.0 or nb == 0.0:
        return False
    return abs(ax * bx + ay * by) <= tol * na * nb


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
    
    def __repr__(self):
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
                EqArea(p1, p2, p3, p4, p5, p6, parents=[self]),
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
                EqArea(p1, p2, p3, p4, p5, p6, parents=[self]),
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

class EqArea(RelationNode):
    def __init__(self, p1: Point, p2: Point, p3: Point, p4: Point, p5: Point, p6: Point,
                 index: Optional[int] = None,
                 parents: Optional[set[RelationNode]] = None,
                 rule: Optional[str] = None):
        super().__init__(
            name="eqarea",
            index=index,
            relation=(frozenset({p1, p2, p3}), frozenset({p4, p5, p6})),
            representation=f"eqarea {p1} {p2} {p3} {p4} {p5} {p6}",
            parents=parents,
            rule=rule,
            equivalent=[]
        )

        self.points = [p1, p2, p3, p4, p5, p6]
