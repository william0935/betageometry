"""Parsing and applying proposed constructions.

Model output is untrusted text, so the parser is the trust boundary between whatever
Gemma emits and the geometry engine. These tests pin down what it must reject.
"""

import pytest

from constructions import Canva
from rabbits import (
    CONSTRUCTION_ARITY,
    ConstructionError,
    RandomProposer,
    apply_construction,
    format_construction,
    parse_construction,
)
from relations import Point

POINTS = [Point("A", 0.0, 0.0), Point("B", 4.0, 0.0), Point("C", 1.0, 3.0)]


def test_parses_a_well_formed_call():
    name, points = parse_construction("foot(A, B, C)", POINTS)
    assert name == "foot"
    assert [p.name for p in points] == ["A", "B", "C"]


def test_tolerates_surrounding_text_and_spacing():
    name, points = parse_construction("  midpoint( A ,B )\n", POINTS)
    assert name == "midpoint"
    assert [p.name for p in points] == ["A", "B"]


@pytest.mark.parametrize("text, reason", [
    ("nonsense", "not a call at all"),
    ("teleport(A, B)", "construction not in the vocabulary"),
    ("foot(A, B)", "wrong arity"),
    ("foot(A, B, C, A)", "wrong arity"),
    ("foot(A, B, Z)", "unknown point"),
    ("foot(A, A, B)", "same point twice"),
])
def test_rejects_bad_output(text, reason):
    with pytest.raises(ConstructionError):
        parse_construction(text, POINTS)


def test_cannot_reach_arbitrary_canva_attributes():
    """The whitelist stops a proposal naming a method that is not a construction."""
    for attr in ("plot", "add_point", "free", "__init__"):
        with pytest.raises(ConstructionError):
            parse_construction(f"{attr}(A, B)", POINTS)


def test_round_trips_through_format():
    text = format_construction("foot", POINTS)
    assert text == "foot(A, B, C)"
    name, points = parse_construction(text, POINTS)
    assert name == "foot" and points == POINTS


def test_apply_returns_nothing_for_a_degenerate_configuration():
    canva = Canva([], {}, {}, {})
    collinear = [Point("A", 0.0, 0.0), Point("B", 1.0, 0.0), Point("C", 2.0, 0.0)]
    points, relations = apply_construction(canva, "foot", collinear, collinear)
    assert points == [] and relations == []


def test_apply_rejects_a_duplicate_of_an_existing_point():
    canva = Canva([], {}, {}, {})
    a, b = Point("A", 0.0, 0.0), Point("B", 4.0, 0.0)
    midpoint = Point("M", 2.0, 0.0)
    points, relations = apply_construction(canva, "midpoint", [a, b], [a, b, midpoint])
    assert points == [], "a construction landing on an existing point adds nothing"


def test_random_proposer_only_emits_parseable_calls():
    proposer = RandomProposer()
    calls = proposer.propose("", POINTS, k=40)
    assert calls
    for call in calls:
        name, points = parse_construction(call, POINTS)
        assert len(points) == CONSTRUCTION_ARITY[name]


def test_random_proposer_weights_must_be_a_distribution():
    with pytest.raises(ValueError):
        RandomProposer(weights={"foot": 0.5, "midpoint": 0.2})
