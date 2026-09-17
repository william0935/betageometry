"""Synthetic training-data generation.

The shape of the task: start from a random configuration of free points, add an
auxiliary point, and ask the engine what it can now prove that it could not prove
before. Any fact that (a) does not mention the auxiliary point but (b) needed it in the
derivation is a case where drawing that point was the key step -- exactly the decision a
model should learn. Each such fact becomes one example:

    input_text : "<premises joined by '; '>; ? <goal>"
    output_text: "foot(A, B, C)"

The premises are the roots of the fact's derivation with the auxiliary point's own
construction relations removed, so the model is asked to invent the point rather than
being handed it.
"""

import json
import math
import os
import random
from itertools import combinations
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from constructions import Canva
from dd_ar import DDWithAR
from problem import Problem
from rabbits import (
    CONSTRUCTION_ARITY,
    RandomProposer,
    apply_construction,
    format_construction,
    parse_construction,
)
from relations import (
    Congruent,
    EqualAngle,
    Midpoint,
    Parallel,
    Perpendicular,
    Point,
    RelationNode,
    angles_close_mod_pi,
    lines_look_parallel,
    lines_look_perpendicular,
)

# Past this many points the deductive closure grows faster than it produces anything
# useful (each simtri fans out into six equivalent relations), so generation stops.
MAX_POINTS = 16


def collect_roots(relation: RelationNode) -> List[RelationNode]:
    """Return the premises a relation ultimately rests on.

    Iterative with a visited set: the derivation graph is a DAG whose nodes are shared
    heavily (a single Cyclic fans out into four `equivalent` children), so a plain
    recursive descent re-walks each node once per distinct path and blows up
    exponentially -- and hits the recursion limit on deep chains.
    """
    roots: List[RelationNode] = []
    visited = set()
    stack = [relation]
    while stack:
        node = stack.pop()
        if id(node) in visited:
            continue
        visited.add(id(node))
        # `parents` defaults to an empty set, so testing it against [] or None never
        # identifies a root and the premise list comes back empty.
        if not node.parents:
            roots.append(node)
            continue
        stack.extend(node.parents)
    return roots


def numeric_goals(points: Sequence[Point], new_points: Sequence[Point]) -> List[RelationNode]:
    """Every relation that holds numerically in the diagram and involves a new point.

    These are handed to the solver as goals so that `check_relation` will try to prove
    them through the AR tables; whatever it proves becomes candidate training data.
    """
    goals: List[RelationNode] = []
    new_ids = {id(p) for p in new_points}

    def involves(*ps):
        return any(id(p) in new_ids for p in ps)

    for p1, p2 in combinations(points, 2):
        for p3, p4 in combinations(points, 2):
            if {id(p1), id(p2)} == {id(p3), id(p4)} or not involves(p1, p2, p3, p4):
                continue
            d1 = math.hypot(p2.x - p1.x, p2.y - p1.y)
            d2 = math.hypot(p4.x - p3.x, p4.y - p3.y)
            if d1 > 0 and math.isclose(d1, d2, rel_tol=1e-6):
                goals.append(Congruent(p1, p2, p3, p4))
            # Scale-invariant cross/dot products. The previous gradient comparisons
            # needed a special case for vertical lines and used
            # `math.isclose(grad, 0.0, rel_tol=...)`, which is only ever true for exactly
            # 0.0 -- rel_tol scales by the magnitude of the operands -- so horizontal and
            # vertical perpendicular pairs were silently never proposed.
            # Only between segments with no endpoint in common. Two segments sharing a
            # point are parallel exactly when the three points are collinear, so such a
            # `para` is a collinearity statement wearing the wrong predicate -- and its
            # `equivalent` relations expand into zero-width angles like
            # "eqangle X1 X2 X1 ...", which pollute the database and the training set.
            if len({id(p1), id(p2), id(p3), id(p4)}) == 4:
                if lines_look_parallel(p1, p2, p3, p4):
                    goals.append(Parallel(p1, p2, p3, p4))
                if lines_look_perpendicular(p1, p2, p3, p4):
                    goals.append(Perpendicular(p1, p2, p3, p4))

    for p1, p2, p3 in combinations(points, 3):
        if not involves(p1, p2, p3):
            continue
        # midp p1 p2 p3 means p1 is the midpoint of p2p3
        mx, my = (p2.x + p3.x) / 2, (p2.y + p3.y) / 2
        if math.isclose(mx, p1.x, abs_tol=1e-9) and math.isclose(my, p1.y, abs_tol=1e-9):
            goals.append(Midpoint(p1, p2, p3))

    # Equal angles, compared the way the engine compares them: directed, modulo pi.
    # Using undirected degrees here proposes goals the engine cannot express.
    angles = []
    for a, b, c in combinations(points, 3):
        for vertex, u, w in ((b, a, c), (a, b, c), (c, a, b)):
            theta = math.atan2(w.y - vertex.y, w.x - vertex.x) - \
                    math.atan2(u.y - vertex.y, u.x - vertex.x)
            theta %= math.pi
            if min(theta, math.pi - theta) < 1e-6:
                continue  # degenerate: the three points are collinear
            angles.append((u, vertex, w, theta))

    for i in range(len(angles)):
        u1, v1, w1, t1 = angles[i]
        for j in range(i + 1, len(angles)):
            u2, v2, w2, t2 = angles[j]
            if not involves(u1, v1, w1, u2, v2, w2):
                continue
            if angles_close_mod_pi(t1, t2, 1e-6):
                goals.append(EqualAngle(u1, v1, w1, u2, v2, w2))

    return goals


def generate_for_problem(problem: Problem,
                         canva: Canva,
                         rounds: int = 4,
                         proposer=None,
                         deduction_iterations: int = 50,
                         verbose: bool = False) -> List[Dict[str, str]]:
    """Add `rounds` auxiliary points, mining training examples after each one."""
    proposer = proposer if proposer is not None else RandomProposer()
    solver = DDWithAR(problem)
    solver.apply_deduction_rules(deduction_iterations, canva, stop_when_solved=False)

    known = {id(r) for rels in problem.relations.values() for r in rels}
    examples: List[Dict[str, str]] = []
    seen_examples = set()

    for round_index in range(rounds):
        if len(problem.points) >= MAX_POINTS:
            if verbose:
                print(f"reached {MAX_POINTS} points; stopping")
            break

        applied = _add_one_rabbit(solver, canva, proposer)
        if applied is None:
            if verbose:
                print("no usable construction this round")
            break
        call, new_points, construction_relations = applied
        if verbose:
            print(f"round {round_index}: {call} -> {[p.name for p in new_points]}")

        for goal in numeric_goals(problem.points, new_points):
            problem.add_goal(goal)

        solver.apply_deduction_rules(deduction_iterations, canva, stop_when_solved=False)

        construction_ids = {id(r) for r in construction_relations}
        new_ids = {id(p) for p in new_points}
        for rels in problem.relations.values():
            for relation in rels:
                if id(relation) in known:
                    continue
                # The interesting facts are the ones the auxiliary point does not appear
                # in -- those are the ones where drawing it actually bought something.
                if any(id(p) in new_ids for p in relation.points):
                    continue
                premises = [r for r in collect_roots(relation)
                            if id(r) not in construction_ids
                            and not any(id(p) in new_ids for p in r.points)]
                if not premises:
                    continue

                input_text = "; ".join(sorted(r.representation for r in premises))
                input_text += f"; ? {relation.representation}"
                key = (input_text, call)
                if key in seen_examples:
                    continue
                seen_examples.add(key)
                examples.append({"input_text": input_text, "output_text": call})

        known = {id(r) for rels in problem.relations.values() for r in rels}

    return examples


def _add_one_rabbit(solver: DDWithAR, canva: Canva, proposer, attempts: int = 25
                    ) -> Optional[Tuple[str, List[Point], List[RelationNode]]]:
    """Ask the proposer for a construction until one applies cleanly."""
    problem = solver.problem
    for _ in range(attempts):
        candidates = proposer.propose("", problem.points, k=1)
        if not candidates:
            continue
        try:
            name, points = parse_construction(candidates[0], problem.points)
        except Exception:
            continue

        new_points, relations = apply_construction(canva, name, points, problem.points)
        if not new_points:
            continue
        if len(problem.points) + len(new_points) > MAX_POINTS:
            continue

        for point in new_points:
            solver.add_constructed_point(point)
        for relation in relations:
            solver.add_constructed_relation(relation)
        return format_construction(name, points), new_points, relations
    return None


def generate_for_seed(seed: int,
                      rounds: int = 4,
                      num_free_points: int = 4,
                      verbose: bool = False) -> List[Dict[str, str]]:
    """Generate examples from one random configuration, reproducibly."""
    rng = random.Random(seed)
    random.seed(seed)
    np.random.seed(seed % (2 ** 32))

    canva = Canva([], {}, {}, {})
    points = [canva.free() for _ in range(num_free_points)]
    # No assumptions and no reachable goal: the engine is being used to enumerate what is
    # true of the configuration, not to prove a particular statement.
    problem = Problem(f"seed_{seed}", points, [], [])
    return generate_for_problem(problem, canva, rounds=rounds,
                                proposer=RandomProposer(rng=rng), verbose=verbose)


def write_examples(examples: List[Dict[str, str]], path: str):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(examples, f, indent=2)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--rounds", type=int, default=4)
    parser.add_argument("--out", default=None, help="output JSON path")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    data = generate_for_seed(args.seed, rounds=args.rounds, verbose=not args.quiet)
    out = args.out or f"data/training_data_{args.seed}.json"
    write_examples(data, out)
    print(f"seed {args.seed}: {len(data)} examples -> {out}")
