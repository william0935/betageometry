"""Constructions must emit only relations that hold in the figure they just drew.

A construction that computes the wrong point while asserting the right-sounding relation
poisons everything downstream: the solver treats the assertion as a premise, so anything
derived from it is unsound, and in data generation it becomes a training target. These
tests fuzz every construction over random configurations and check each emitted relation
against the coordinates.
"""

import math
import random

import pytest

from constructions import Canva
from geometry_check import check
from rabbits import CONSTRUCTION_ARITY
from relations import Point

CONSTRUCTIONS = sorted(CONSTRUCTION_ARITY)


def random_points(rng, n=3, span=6.0):
    return [Point(chr(ord("A") + i), rng.uniform(-span, span), rng.uniform(-span, span))
            for i in range(n)]


@pytest.mark.parametrize("name", CONSTRUCTIONS)
def test_construction_emits_only_true_relations(name):
    rng = random.Random(hash(name) & 0xFFFF)
    canva = Canva([], {}, {}, {})
    checked = 0
    for _ in range(60):
        points = random_points(rng, CONSTRUCTION_ARITY[name])
        result, relations = getattr(canva, name)(*points)
        if not result:
            continue
        for relation in relations:
            verdict = check(relation)
            if verdict is None:
                continue
            checked += 1
            assert verdict, (
                f"{name}{tuple(p.name for p in points)} emitted a false relation: "
                f"{relation.representation}"
            )
    assert checked > 0, f"{name} never produced a checkable relation"


@pytest.mark.parametrize("name", CONSTRUCTIONS)
def test_construction_never_raises(name):
    """Degenerate inputs must return an empty result, not blow up."""
    rng = random.Random(1234 + len(name))
    canva = Canva([], {}, {}, {})
    for _ in range(40):
        points = random_points(rng, CONSTRUCTION_ARITY[name])
        getattr(canva, name)(*points)  # must not raise


def test_orthocenter_is_the_orthocenter():
    canva = Canva([], {}, {}, {})
    a, b, c = Point("A", 0.0, 0.0), Point("B", 3.0, 4.0), Point("C", 6.0, 0.0)
    h, _ = canva.orthocenter(a, b, c)
    assert h is not None
    assert math.isclose(h.x, 3.0, abs_tol=1e-9)
    assert math.isclose(h.y, 2.25, abs_tol=1e-9)


def test_intersect_lines_handles_vertical_and_parallel():
    canva = Canva([], {}, {}, {})
    # a vertical line: the slope form used to divide by zero here
    p, _ = canva.intersect_lines(Point("A", 0.0, 0.0), Point("B", 0.0, 4.0),
                                 Point("C", 2.0, 1.0), Point("D", 5.0, 1.0))
    assert p is not None
    assert math.isclose(p.x, 0.0, abs_tol=1e-9) and math.isclose(p.y, 1.0, abs_tol=1e-9)

    # genuinely parallel lines have no intersection to report
    q, relations = canva.intersect_lines(Point("P", 0.0, 0.0), Point("Q", 1.0, 1.0),
                                         Point("R", 0.0, 1.0), Point("S", 1.0, 2.0))
    assert q is None and relations == []


def test_tangent_points_lie_on_the_circle():
    canva = Canva([], {}, {}, {})
    a, o, b = Point("A", 10.0, 0.0), Point("O", 0.0, 0.0), Point("B", 3.0, 0.0)
    (t1, t2), _ = canva.tangent2(a, o, b)
    radius = math.hypot(b.x - o.x, b.y - o.y)
    for t in (t1, t2):
        assert math.isclose(math.hypot(t.x - o.x, t.y - o.y), radius, rel_tol=1e-9)
        # tangency: the radius to the point is perpendicular to the line from a
        dot = (t.x - o.x) * (t.x - a.x) + (t.y - o.y) * (t.y - a.y)
        assert abs(dot) < 1e-6


def test_auxiliary_names_avoid_collisions_with_the_diagram():
    canva = Canva([], {"X1": (0.0, 0.0)}, {}, {})
    created = canva.add_point(1.0, 1.0)
    assert created.name != "X1"
