"""Inequality reasoning.

The central property is that inequality rows may only be *added*, never scaled by a
negative number. A span-based test -- which is what the equality tables use -- concludes
`b >= a` from `a >= b`, so most of these tests are about the unsound direction staying
unprovable.
"""

import math

import pytest

from ar_inequalities import (
    AngleInequalityTable,
    InequalityTable,
    SegmentInequalityTable,
    nonneg_combination,
)
from constructions import Canva
from dd_ar import DDWithAR
from problem import Problem
from relations import (
    Congruent,
    GreaterAngle,
    GreaterEqSegment,
    GreaterSegment,
    Point,
    RelationNode,
)

SOURCE = RelationNode("gtseg", 1, ("k", 1))


def make_table(columns, inequalities=(), equalities=()):
    table = InequalityTable(list(columns))
    for row, strict in inequalities:
        table.add_inequality(row, SOURCE, strict=strict)
    for row in equalities:
        table.add_equality(row, SOURCE)
    return table


# --------------------------------------------------------------------------- solver

def test_only_addition_reversal_is_not_derivable():
    """The whole point: scaling an inequality by -1 must not be allowed."""
    table = make_table(["a", "b"], inequalities=[([1, -1], True)])
    assert table.implies([1, -1])[0], "a > b should follow from itself"
    assert not table.implies([-1, 1])[0], "b > a is the opposite claim"


def test_transitivity_holds_but_only_forwards():
    table = make_table(["a", "b", "c"],
                       inequalities=[([1, -1, 0], True), ([0, 1, -1], True)])
    assert table.implies([1, 0, -1])[0], "a > b > c gives a > c"
    assert not table.implies([-1, 0, 1])[0], "c > a does not follow"


def test_opposing_inequalities_are_both_kept():
    """`a >= b` and `b >= a` are linearly dependent but jointly mean a == b.

    The equality tables drop dependent rows as redundant; here that would lose the
    second premise entirely.
    """
    table = make_table(["a", "b"],
                       inequalities=[([1, -1], False), ([-1, 1], False)])
    assert table.implies([1, -1], strict=False)[0]
    assert table.implies([-1, 1], strict=False)[0]


def test_equalities_may_be_used_with_either_sign():
    table = make_table(["a", "b", "c"],
                       inequalities=[([0, 1, -1], True)],
                       equalities=[[1, -1, 0]])
    assert table.implies([1, 0, -1])[0], "a == b and b > c gives a > c"
    assert not table.implies([-1, 0, 1])[0]


def test_strict_conclusion_needs_a_strict_premise():
    table = make_table(["a", "b", "c"],
                       inequalities=[([1, -1, 0], False), ([0, 1, -1], False)])
    assert table.implies([1, 0, -1], strict=False)[0], "a >= b >= c gives a >= c"
    assert not table.implies([1, 0, -1], strict=True)[0], "but not a > c"


def test_empty_table_proves_nothing():
    assert not InequalityTable(["a", "b"]).implies([1, -1])[0]


def test_certificate_names_the_rows_actually_used():
    r1 = RelationNode("gtseg", 1, ("k", 1))
    r2 = RelationNode("gtseg", 2, ("k", 2))
    r3 = RelationNode("gtseg", 3, ("k", 3))
    table = InequalityTable(["a", "b", "c", "d"])
    table.add_inequality([1, -1, 0, 0], r1, strict=True)
    table.add_inequality([0, 1, -1, 0], r2, strict=True)
    table.add_inequality([0, 0, 0, 1], r3, strict=True)   # irrelevant to the query
    proved, used = table.implies([1, 0, -1, 0])
    assert proved
    assert used == {r1, r2}, "the unrelated row must not appear in the proof"


def test_solver_agrees_with_brute_force():
    import itertools
    import random

    import numpy as np

    def brute(nonneg, free, target, lim=3):
        for lam in itertools.product(range(0, lim + 1), repeat=len(nonneg)):
            for mu in itertools.product(range(-lim, lim + 1), repeat=len(free)):
                v = np.zeros(len(target))
                for c, row in zip(lam, nonneg):
                    v += c * np.array(row, dtype=float)
                for c, row in zip(mu, free):
                    v += c * np.array(row, dtype=float)
                if np.allclose(v, target):
                    return True
        return False

    rng = random.Random(0)
    for _ in range(150):
        ncol = rng.randint(2, 4)
        nonneg = [[rng.randint(-2, 2) for _ in range(ncol)]
                  for _ in range(rng.randint(1, 3))]
        free = [[rng.randint(-2, 2) for _ in range(ncol)]
                for _ in range(rng.randint(0, 2))]
        target = [rng.randint(-3, 3) for _ in range(ncol)]
        found = nonneg_combination(nonneg, free, target) is not None
        # The LP may use fractional coefficients the bounded integer search cannot
        # reach, so it can only ever find *more*; it must never find less.
        assert found or not brute(nonneg, free, target), \
            f"solver missed a combination: {nonneg} {free} {target}"


# ------------------------------------------------------------------- lazy activation

def scalene():
    """|AB| = 6, |BC| ~ 5.83, |AC| ~ 3.16 -- all three sides distinct."""
    return Point("A", 0.0, 0.0), Point("B", 6.0, 0.0), Point("C", 1.0, 3.0)


def build(points, assumptions, goals):
    problem = Problem("ineq", list(points), assumptions, goals)
    canva = Canva(list(points), {p.name: (p.x, p.y) for p in points}, {}, {})
    return problem, canva, DDWithAR(problem)


def test_inequality_machinery_is_off_without_an_inequality():
    a, b, c = scalene()
    _, _, solver = build((a, b, c), [Congruent(a, b, a, c)], [])
    assert solver.uses_inequalities is False
    assert solver.segment_ineq_table is None
    assert solver.angle_ineq_table is None
    rule_names = {r.__name__ for r in solver.rules}
    assert not any("gt" in n for n in rule_names), "inequality rules must not be loaded"


@pytest.mark.parametrize("where", ["assumption", "goal"])
def test_an_inequality_anywhere_turns_it_on(where):
    a, b, c = scalene()
    rel = GreaterSegment(a, b, a, c)
    assumptions = [rel] if where == "assumption" else []
    goals = [rel] if where == "goal" else []
    _, _, solver = build((a, b, c), assumptions, goals)
    assert solver.uses_inequalities is True
    assert solver.segment_ineq_table is not None


def test_existing_problems_are_untouched(solved_problems):
    """No stock problem mentions an inequality, so none should switch the mode on."""
    from conftest import load_problem

    for name in ("problem1", "problem3", "usamo_2023_p1"):
        problem, canva = load_problem(name)
        assert DDWithAR(problem).uses_inequalities is False, name


# ------------------------------------------------------------------------ deduction

def solve_ineq(points, assumptions, goals, iterations=20):
    problem, canva, solver = build(points, assumptions, goals)
    solver.apply_deduction_rules(iterations, canva)
    return problem


def test_larger_side_faces_larger_angle():
    a, b, c = scalene()
    problem = solve_ineq((a, b, c), [GreaterSegment(a, b, a, c)],
                         [GreaterAngle(a, c, b, a, b, c)])
    assert problem.solved


def test_larger_angle_faces_larger_side():
    a, b, c = scalene()
    problem = solve_ineq((a, b, c), [GreaterAngle(a, c, b, a, b, c)],
                         [GreaterSegment(a, b, a, c)])
    assert problem.solved


def test_transitivity_across_stated_inequalities():
    a, b, c = scalene()
    d = Point("D", 2.0, 1.0)
    problem = solve_ineq((a, b, c, d),
                         [GreaterSegment(a, b, a, c), GreaterSegment(a, c, a, d)],
                         [GreaterSegment(a, b, a, d)])
    assert problem.solved


def test_congruence_combines_with_an_inequality():
    a = Point("A", 0.0, 0.0)
    b = Point("B", 6.0, 0.0)
    d = Point("D", 3.0, 0.0)     # |AD| = 3
    e = Point("E", 0.0, 3.0)     # |AE| = 3, so AD and AE are congruent
    problem = solve_ineq((a, b, d, e),
                         [Congruent(a, d, a, e), GreaterSegment(a, b, a, d)],
                         [GreaterSegment(a, b, a, e)])
    assert problem.solved
    proof = problem.trace_back()
    assert "cong A D A E" in proof, "the congruence must appear in the proof"


def test_the_reverse_inequality_is_never_proved():
    a, b, c = scalene()
    problem = solve_ineq((a, b, c), [GreaterSegment(a, b, a, c)],
                         [GreaterSegment(a, c, a, b)])
    assert not problem.solved, "proving the reverse would mean negative coefficients"


def test_strict_premise_proves_the_non_strict_goal():
    a, b, c = scalene()
    problem = solve_ineq((a, b, c), [GreaterSegment(a, b, a, c)],
                         [GreaterEqSegment(a, b, a, c)])
    assert problem.solved


def test_non_strict_premise_does_not_prove_the_strict_goal():
    a, b, c = scalene()
    problem = solve_ineq((a, b, c), [GreaterEqSegment(a, b, a, c)],
                         [GreaterSegment(a, b, a, c)])
    assert not problem.solved


def test_nothing_derived_contradicts_the_diagram():
    """Every inequality the engine produces must hold of the actual coordinates."""
    from relations import angle_magnitude

    a, b, c = scalene()
    d = Point("D", 2.0, 1.0)
    problem = solve_ineq((a, b, c, d),
                         [GreaterSegment(a, b, a, c), GreaterAngle(a, c, b, a, b, c)],
                         [])

    def length(p, q):
        return math.hypot(q.x - p.x, q.y - p.y)

    for name in ("gtseg", "gteseg"):
        for rel in problem.relations[name]:
            p1, p2, p3, p4 = rel.points
            lhs, rhs = length(p1, p2), length(p3, p4)
            assert lhs > rhs - 1e-9, f"{rel.representation} is false in the diagram"

    for name in ("gtangle", "gteangle"):
        for rel in problem.relations[name]:
            p1, p2, p3, p4, p5, p6 = rel.points
            lhs = angle_magnitude(p1, p2, p3)
            rhs = angle_magnitude(p4, p5, p6)
            assert lhs > rhs - 1e-9, f"{rel.representation} is false in the diagram"


def test_comparing_a_quantity_with_itself_is_rejected():
    a, b, c = scalene()
    problem, _, _ = build((a, b, c), [], [])
    assert problem.add_relation(GreaterSegment(a, b, a, b)) is None, \
        "|AB| > |AB| is false and must never be stored"
    assert problem.add_relation(GreaterEqSegment(a, b, a, b)) is None, \
        "|AB| >= |AB| is vacuous"


def test_statement_dsl_round_trip():
    from read_in_relations import get_relation

    points = [Point(n, 0.0, 0.0) for n in "ABCDEF"]
    for text, expected in (
        ("gtseg A B C D", "gtseg"),
        ("gteseg A B C D", "gteseg"),
        ("gtangle A B C D E F", "gtangle"),
        ("gteangle A B C D E F", "gteangle"),
    ):
        relation = get_relation(text, points)
        assert relation.name == expected
        assert relation.representation == text
