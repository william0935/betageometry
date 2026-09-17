"""Solve-with-auxiliary-points loop.

`DDWithAR.apply_deduction_rules` closes everything derivable from the points already on
the diagram. When that is not enough, the loop here asks a `RabbitProposer` for an extra
point, adds it, and runs the engine again -- the neural half of the AlphaGeometry-style
split, with the symbolic half unchanged and still responsible for every proof step.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Sequence

from constructions import Canva
from dd_ar import DDWithAR
from problem import Problem
from rabbits import (
    ConstructionError,
    RabbitProposer,
    apply_construction,
    parse_construction,
)
from relations import RelationNode


def problem_to_text(problem: Problem) -> str:
    """Render a problem in the DSL the proposer is trained on.

    Matches the `text_files/*.txt` statement format and the `input_text` field of the
    generated training data: assumptions separated by "; ", then "? " and the goal.
    """
    parts = [r.representation for r in problem.assumptions]
    goals = " ".join(f"? {g.representation}" for g in problem.remaining_goals)
    body = "; ".join(parts)
    return f"{body}; {goals}" if body else goals


@dataclass
class SearchResult:
    solved: bool
    constructions: List[str] = field(default_factory=list)
    rounds: int = 0


def solve(problem: Problem,
          canva: Canva,
          proposer: Optional[RabbitProposer] = None,
          max_rounds: int = 8,
          candidates_per_round: int = 4,
          deduction_iterations: int = 50,
          solver: Optional[DDWithAR] = None,
          verbose: bool = True) -> SearchResult:
    """Run the engine, adding proposed auxiliary points until the goal closes.

    With `proposer=None` this is exactly the plain deductive solve, so the symbolic path
    stays available and testable on its own.
    """
    solver = solver if solver is not None else DDWithAR(problem)

    if solver.apply_deduction_rules(deduction_iterations, canva):
        return SearchResult(solved=True, rounds=0)

    if proposer is None:
        return SearchResult(solved=False, rounds=0)

    applied: List[str] = []
    for round_index in range(max_rounds):
        text = problem_to_text(problem)
        try:
            candidates = proposer.propose(text, problem.points, k=candidates_per_round)
        except Exception as exc:  # a proposer failure must not lose the partial proof
            if verbose:
                print(f"[search] proposer failed: {type(exc).__name__}: {exc}")
            break

        progressed = False
        for candidate in candidates:
            try:
                name, points = parse_construction(candidate, problem.points)
            except ConstructionError as exc:
                if verbose:
                    print(f"[search] rejected {candidate!r}: {exc}")
                continue

            new_points, relations = apply_construction(
                canva, name, points, problem.points
            )
            if not new_points:
                if verbose:
                    print(f"[search] {candidate} was degenerate here, skipping")
                continue

            for point in new_points:
                solver.add_constructed_point(point)
            for relation in relations:
                solver.add_constructed_relation(relation)

            applied.append(candidate)
            progressed = True
            if verbose:
                print(f"[search] round {round_index}: applied {candidate} "
                      f"-> {[p.name for p in new_points]}")
            break

        if not progressed:
            if verbose:
                print("[search] no candidate could be applied; stopping")
            break

        if solver.apply_deduction_rules(deduction_iterations, canva):
            return SearchResult(solved=True, constructions=applied,
                                rounds=round_index + 1)

    return SearchResult(solved=problem.is_solved(), constructions=applied,
                        rounds=len(applied))
