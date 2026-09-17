"""End-to-end checks on the deductive engine.

Two independent properties:

* **Completeness** -- the problems that were solvable stay solvable. A rule change that
  loses a deduction shows up as a named failure.
* **Soundness** -- every relation the engine derives is actually true of the diagram it
  was derived from. This is the check that catches transposed indices and miscomputed
  constructions, which the symbolic layer cannot detect on its own.
"""

import pytest

from conftest import ALL_PROBLEMS, SOLVABLE, solve_problem
from geometry_check import check, false_relations


@pytest.mark.parametrize("name", sorted(SOLVABLE))
def test_solvable_problems_still_solve(name, solved_problems):
    assert solved_problems[name].solved, f"{name} no longer closes"


@pytest.mark.parametrize("name", sorted(ALL_PROBLEMS))
def test_derived_relations_are_true_in_the_diagram(name, solved_problems):
    problem = solved_problems[name]
    bad = false_relations(problem)
    assert not bad, (
        f"{name} derived {len(bad)} relation(s) contradicted by the diagram: "
        + "; ".join(f"{r.representation} (via {r.rule or 'assumption'})" for r in bad[:5])
    )


@pytest.mark.parametrize("name", sorted(ALL_PROBLEMS))
def test_every_relation_is_checked_or_knowingly_unchecked(name, solved_problems):
    """Guards the guard: if a relation type stops being checkable the suite says so."""
    problem = solved_problems[name]
    unchecked = {r.name for rels in problem.relations.values() for r in rels
                 if check(r) is None}
    assert not unchecked - {"eqratio", "simtri1", "simtri2", "cyclic"}, (
        f"{name} has relations no numeric check covers: {unchecked}"
    )


def test_goal_is_actually_derived_and_the_proof_is_non_empty():
    problem = solve_problem("problem1")
    assert problem.solved
    # The proof may cite an equivalent ordering of the goal ("contri2 B R Y D Y R" for
    # "contri2 R B D Y D B"), so check the canonical relation rather than the wording.
    for goal in problem.goals:
        assert goal.relation in problem.relation_keys[goal.name], \
            f"{goal.representation} was never derived"
    assert problem.trace_back().strip(), "solved problem produced an empty proof"


def test_repr_is_stable_across_calls():
    """__repr__ runs trace_back, which used to renumber relations in place."""
    problem = solve_problem("problem1")
    assert repr(problem) == repr(problem)


def test_solver_without_proposer_matches_plain_deduction():
    from search import solve
    from conftest import load_problem

    problem, canva = load_problem("problem1")
    result = solve(problem, canva, proposer=None, verbose=False)
    assert result.solved
    assert result.constructions == []
