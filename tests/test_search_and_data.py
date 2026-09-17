"""The LLM-in-the-loop search, and the data generator that trains it.

The search is exercised with a scripted proposer rather than a real model, so these run
offline and deterministically. What is being tested is the plumbing around the model --
whether a proposed construction reaches the solver, whether bad output is survivable --
not the model's taste in constructions.
"""

import json

import pytest

from conftest import load_problem
from gemma import INSTRUCTION, SEPARATOR, build_example, build_prompt
from geometry_check import check
from search import problem_to_text, solve


class ScriptedProposer:
    """Returns a fixed list of candidates, recording what it was asked."""

    def __init__(self, *responses):
        self.responses = list(responses)
        self.prompts = []

    def propose(self, problem_text, points, k=1):
        self.prompts.append(problem_text)
        return self.responses.pop(0) if self.responses else []


def test_search_applies_a_proposed_construction():
    problem, canva = load_problem("test2")  # not solvable by deduction alone
    before = len(problem.points)
    names = [p.name for p in problem.points]
    proposer = ScriptedProposer([f"midpoint({names[0]}, {names[1]})"])

    solve(problem, canva, proposer=proposer, max_rounds=1, verbose=False)

    assert len(problem.points) == before + 1, "the new point never reached the solver"
    assert proposer.prompts, "the proposer was never consulted"


def test_search_survives_unparseable_output():
    problem, canva = load_problem("test2")
    names = [p.name for p in problem.points]
    proposer = ScriptedProposer(
        ["not a construction", "teleport(A, B)", f"midpoint({names[0]}, {names[1]})"]
    )
    result = solve(problem, canva, proposer=proposer, max_rounds=1, verbose=False)
    assert result.constructions == [f"midpoint({names[0]}, {names[1]})"]


def test_search_survives_a_proposer_that_raises():
    class Broken:
        def propose(self, *a, **k):
            raise RuntimeError("model died")

    problem, canva = load_problem("test2")
    result = solve(problem, canva, proposer=Broken(), max_rounds=2, verbose=False)
    assert result.solved is False  # degrades to the plain deductive answer


def test_search_without_a_proposer_is_plain_deduction():
    problem, canva = load_problem("problem1")
    result = solve(problem, canva, proposer=None, verbose=False)
    assert result.solved and result.constructions == []


def test_problem_text_matches_the_statement_dsl():
    problem, _ = load_problem("problem1")
    text = problem_to_text(problem)
    assert "? " in text, "the goal must be marked with '?'"
    for assumption in problem.assumptions:
        assert assumption.representation in text


def test_prompt_is_shared_between_training_and_inference():
    """Training and inference must build byte-identical prompts."""
    statement = "cong A B A C; ? perp A D B C"
    prompt = build_prompt(statement)
    example = build_example(statement, "foot(A, B, C)", "<eos>")
    assert example.startswith(prompt), "training text and inference prompt have diverged"
    assert example[len(prompt):] == "foot(A, B, C)<eos>"
    assert INSTRUCTION in prompt and prompt.endswith(SEPARATOR)


def test_generated_examples_are_well_formed(tmp_path):
    import data_generation as dg

    examples = []
    for seed in range(12):
        examples.extend(dg.generate_for_seed(seed, rounds=3))
    assert examples, "generator produced nothing across 12 seeds"

    from rabbits import CONSTRUCTION_ARITY, parse_construction
    from relations import Point

    for example in examples:
        assert set(example) == {"input_text", "output_text"}
        assert "; ? " in example["input_text"], example["input_text"]
        # the target must be a construction the solver can actually execute
        name = example["output_text"].split("(")[0]
        assert name in CONSTRUCTION_ARITY, example["output_text"]

    path = tmp_path / "training_data.json"
    dg.write_examples(examples, str(path))
    assert len(json.loads(path.read_text())) == len(examples)


def test_generated_goals_are_true_in_their_diagram():
    """A training goal that is false teaches the model to chase impossible statements."""
    import data_generation as dg
    from constructions import Canva
    from problem import Problem
    from rabbits import RandomProposer
    import random

    for seed in range(8):
        rng = random.Random(seed)
        random.seed(seed)
        canva = Canva([], {}, {}, {})
        points = [canva.free() for _ in range(4)]
        problem = Problem(f"seed_{seed}", points, [], [])
        dg.generate_for_problem(problem, canva, rounds=3,
                                proposer=RandomProposer(rng=rng))
        bad = [r for rels in problem.relations.values() for r in rels
               if check(r) is False]
        assert not bad, (
            f"seed {seed} derived false relations: "
            + "; ".join(r.representation for r in bad[:5])
        )
