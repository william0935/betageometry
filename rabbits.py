"""Auxiliary-point ("rabbit") proposal.

The deductive engine in `dd_ar.py` is complete for what it can see, but many problems
only close once an extra point is drawn. Choosing that point is the part that is not
mechanical, so it is the part a language model is asked to do.

This module is the seam between the two halves:

  * `CONSTRUCTIONS` is the vocabulary -- the construction calls a proposer may emit, and
    the arity of each.
  * `parse_construction` turns the text a model produces ("foot(A, B, C)") back into a
    call against a `Canva`.
  * `RabbitProposer` is the interface; `RandomProposer` samples from the vocabulary and
    `GemmaProposer` asks a fine-tuned Gemma. Data generation uses the former, solving
    uses the latter, and both go through `apply_construction` so a proposed point
    reaches the solver the same way regardless of where it came from.
"""

import math
import random
import re
from typing import Callable, Dict, List, Optional, Protocol, Sequence, Tuple

from constructions import Canva
from relations import Point, RelationNode


# name -> number of existing points the construction consumes.
# This is both the sampling vocabulary and the whitelist used when parsing model output:
# a proposer can only ever name a construction that appears here.
CONSTRUCTION_ARITY: Dict[str, int] = {
    "midpoint": 2,
    "mirror": 2,
    "tangent": 2,
    "on_dia": 2,
    "anglebisector": 3,
    "foot": 3,
    "circle": 3,
    "incenter": 3,
    "incenter2": 3,
    "excenter": 3,
    "excenter2": 3,
    "centroid": 3,
    "orthocenter": 3,
    "orthocenter2": 3,
    "reflect": 3,
    "parallel": 3,
    "tangent2": 3,
}

# Sampling weights for RandomProposer. Cheap, frequently-useful constructions are
# weighted up; the compound ones (incenter2 and friends) introduce seven points at once,
# which inflates the diagram fast, so they stay rare.
DEFAULT_WEIGHTS: Dict[str, float] = {
    "midpoint": 0.15,
    "anglebisector": 0.15,
    "foot": 0.15,
    "circle": 0.15,
    "incenter": 0.05,
    "excenter": 0.05,
    "centroid": 0.05,
    "orthocenter": 0.05,
    "reflect": 0.05,
    "parallel": 0.05,
    "tangent": 0.05,
    "tangent2": 0.05,
}

_CALL_RE = re.compile(r"([A-Za-z_][A-Za-z_0-9]*)\s*\(([^)]*)\)")


class ConstructionError(ValueError):
    """Raised when a proposed construction cannot be turned into a real call."""


def parse_construction(text: str, points: Sequence[Point]) -> Tuple[str, List[Point]]:
    """Parse ``name(P, Q, ...)`` into a construction name and the points it names.

    Model output is untrusted text: it may name a construction that does not exist, pass
    the wrong number of arguments, or reference a point that is not in the diagram. Each
    of those raises `ConstructionError` rather than being coerced into something that
    would silently produce a wrong diagram. Only names in `CONSTRUCTION_ARITY` are
    accepted, so a parse can never reach arbitrary attributes on `Canva`.
    """
    match = _CALL_RE.search(text.strip())
    if not match:
        raise ConstructionError(f"no construction call found in {text!r}")

    name = match.group(1)
    if name not in CONSTRUCTION_ARITY:
        raise ConstructionError(f"unknown construction {name!r}")

    raw_args = [a.strip() for a in match.group(2).split(",") if a.strip()]
    expected = CONSTRUCTION_ARITY[name]
    if len(raw_args) != expected:
        raise ConstructionError(
            f"{name} takes {expected} points, got {len(raw_args)}: {raw_args}"
        )

    by_name = {p.name: p for p in points}
    resolved = []
    for arg in raw_args:
        if arg not in by_name:
            raise ConstructionError(f"{name} references unknown point {arg!r}")
        resolved.append(by_name[arg])

    if len({id(p) for p in resolved}) != len(resolved):
        raise ConstructionError(f"{name} names the same point twice: {raw_args}")

    return name, resolved


def format_construction(name: str, points: Sequence[Point]) -> str:
    """Render a construction call in the form used for training targets."""
    return f"{name}({', '.join(p.name for p in points)})"


def apply_construction(
    canva: Canva,
    name: str,
    points: Sequence[Point],
    existing_points: Sequence[Point],
) -> Tuple[List[Point], List[RelationNode]]:
    """Run a construction and return the points and relations it produced.

    Returns ``([], [])`` when the construction is degenerate for these inputs (the
    constructions signal that as either ``(None, [])`` or ``([], [])``) or when it lands
    on a point the diagram already has -- a duplicate point adds no information and gives
    the solver a second name for something it can already talk about.
    """
    func = getattr(canva, name)
    new_points, relations = func(*points)

    if not new_points:
        return [], []
    if not isinstance(new_points, list):
        new_points = [new_points]

    for new_point in new_points:
        for old_point in existing_points:
            if math.isclose(new_point.x, old_point.x, abs_tol=1e-9) and \
               math.isclose(new_point.y, old_point.y, abs_tol=1e-9):
                return [], []

    return new_points, relations


class RabbitProposer(Protocol):
    """Proposes auxiliary constructions for a problem."""

    def propose(self, problem_text: str, points: Sequence[Point], k: int = 1) -> List[str]:
        """Return up to `k` construction calls, best first."""
        ...


class RandomProposer:
    """Samples uniformly-at-random constructions from the weighted vocabulary.

    This is the proposer used for data generation, and the baseline the fine-tuned model
    is measured against.
    """

    def __init__(self, weights: Optional[Dict[str, float]] = None,
                 rng: Optional[random.Random] = None):
        self.weights = dict(weights if weights is not None else DEFAULT_WEIGHTS)
        total = sum(self.weights.values())
        if not math.isclose(total, 1.0, rel_tol=1e-9):
            raise ValueError(f"weights must sum to 1.0 (got {total})")
        self.rng = rng if rng is not None else random

    def propose(self, problem_text: str, points: Sequence[Point], k: int = 1) -> List[str]:
        names = list(self.weights)
        weights = [self.weights[n] for n in names]
        out = []
        for _ in range(k):
            # Re-draw from the full distribution each time. Falling through to the next
            # entry on a failure would hand the failed construction's probability mass to
            # whichever one happened to be listed after it.
            name = self.rng.choices(names, weights=weights, k=1)[0]
            arity = CONSTRUCTION_ARITY[name]
            if len(points) < arity:
                continue
            chosen = self.rng.sample(list(points), arity)
            out.append(format_construction(name, chosen))
        return out
