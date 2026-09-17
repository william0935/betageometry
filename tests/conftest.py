import os
import sys

import matplotlib
matplotlib.use("Agg")  # tests must not try to open a window

# The project is a flat set of top-level modules (the Oscar workflow uploads them into a
# single directory), so the repository root has to be importable from tests/.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest

ALL_PROBLEMS = [
    "problem1", "problem2", "problem3", "problem4", "problem5", "problem6",
    "problem7", "problem8", "problem9", "problem10", "problem11",
    "test", "test2", "test_SAS",
    "usamo_2023_p1", "Yasinsky_Geometry_Olympiad_2023_VIII_p1",
]

# Problems the deductive engine closes without auxiliary points. Locked in so a
# regression shows up as a named failure rather than a silently smaller number.
SOLVABLE = {
    "problem1", "problem2", "problem3", "problem4", "problem5", "problem6",
    "problem7", "problem8", "problem9", "problem10", "problem11",
    "test", "test_SAS",
}


def load_problem(name):
    """Build (problem, canva) for a named problem, without solving it."""
    from constructions import Canva
    from problem import Problem
    from read_in_geogebra_file import parse_picture
    from read_in_relations import read_in_relations
    from relations import Point

    points_dict, lines, circles = parse_picture(f"{name}.ggb")
    points = [Point(n, x, y) for n, (x, y) in points_dict.items()]
    canva = Canva(points, points_dict, lines, circles)
    assumptions, goals = read_in_relations(f"{name}.txt", points)
    return Problem(name, points, assumptions, goals), canva


def solve_problem(name, iterations=50):
    """Run the engine to completion and return the solved Problem."""
    from dd_ar import DDWithAR

    problem, canva = load_problem(name)
    DDWithAR(problem).apply_deduction_rules(iterations, canva)
    return problem


@pytest.fixture(scope="session")
def solved_problems():
    """Every problem solved once, shared across the tests that inspect the results."""
    return {name: solve_problem(name) for name in ALL_PROBLEMS}
