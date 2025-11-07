import matplotlib
import pytest
import numpy as np
import math
from matplotlib import pyplot as plt

from relations import Point
from Problem import Problem
from dd_ar import DDWithAR
from Constructions import ConstructionPlotter

def plotter_ax():
    plotter = ConstructionPlotter()
    fig, ax = plotter.setup_geometry_plot()
    return plotter, fig, ax

def test_find_intersection_and_counter(plotter_ax):
    plotter, fig, ax = plotter_ax
    A = Point("A", 0.0, 0.0)
    B = Point("B", 1.0, 1.0)
    C = Point("C", 0.0, 1.0)
    D = Point("D", 1.0, 0.0)
    start_counter = plotter.counter
    X = plotter.find_intersection(A, B, C, D)
    assert isinstance(X, Point)
    assert X.name.startswith("X")
    assert abs(X.x - 0.5) < 1e-8 and abs(X.y - 0.5) < 1e-8
    assert plotter.counter == start_counter + 1

def test_plot_midpoint(plotter_ax):
    plotter, fig, ax = plotter_ax
    A = Point("A", 0.0, 0.0)
    B = Point("B", 2.0, 0.0)
    M = plotter.plot_midpoint(ax, "M_AB", A, B)
    assert isinstance(M, Point)
    assert M.name == "M_AB"
    assert abs(M.x - 1.0) < 1e-8 and abs(M.y - 0.0) < 1e-8

def test_plotting_helpers_do_not_raise(plotter_ax):
    plotter, fig, ax = plotter_ax
    # plotting dictionaries
    points_dict = {"P": (0.0, 0.0), "Q": (1.0, 0.0)}
    lines_eq = {"l1": (0.0, 1.0, 0.0)}   # y = -c/b (simple)
    circles = {"c1": (0.0, 0.0, 1.0)}
    # call plot helpers
    plotter.plot_points(ax, points_dict)
    plotter.plot_lines_from_eq(ax, lines_eq)
    plotter.plot_circles_from_eq(ax, circles)

    # new primitives
    A = Point("A", 0.0, 0.0)
    B = Point("B", 1.0, 0.0)
    plotter.plot_newpoint(ax, A)
    plotter.plot_newline(ax, A, B)
    plotter.plot_newcircle(ax, A, B)  # circle centered at A through B

def test_anglebisector_and_tables(plotter_ax):
    plotter, fig, ax = plotter_ax
    A = Point("A", 0.0, 0.0)
    B = Point("B", 1.0, 0.0)
    C = Point("C", 0.0, 1.0)
    prob = Problem("t", [A, B, C], assumptions=[], goals=[])
    dd = DDWithAR(prob)
    new_rels = plotter.plot_anglebisector(ax, dd, A, B, C)
    # should produce EqualAngle and Collinear relations
    assert any(rel.__class__.__name__ in ("EqualAngle", "Collinear") for rel in new_rels)
    # angle_table should have columns for a name like "AX{n}" etc (strings added by the method)
    # check that at least one of the expected string-keys exists
    expected_keys = {A.name + "X" + str(plotter.counter - 1), B.name + "X" + str(plotter.counter - 1)}
    assert any(k in dd.angle_table.col_id for k in expected_keys)

def test_incenter_and_circumcenter(plotter_ax):
    plotter, fig, ax = plotter_ax
    A = Point("A", 0.0, 0.0)
    B = Point("B", 2.0, 0.0)
    C = Point("C", 0.0, 2.0)
    prob = Problem("t", [A, B, C], assumptions=[], goals=[])
    dd = DDWithAR(prob)
    in_rels = plotter.incenter(ax, dd, A, B, C)
    assert isinstance(in_rels, list)
    circ_rels = plotter.circumcenter(ax, dd, A, B, C)
    assert isinstance(circ_rels, list)

def test_plot_perp_updates_tables(plotter_ax):
    plotter, fig, ax = plotter_ax
    A = Point("A", 0.0, 0.0)
    B = Point("B", 1.0, 0.0)
    P = Point("P", 0.0, 1.0)
    prob = Problem("t", [A, B, P], assumptions=[], goals=[])
    dd = DDWithAR(prob)
    new_rels = plotter.plot_perp(ax, dd, A, B, P)
    assert any(rel.__class__.__name__ == "Perpendicular" for rel in new_rels)

if __name__ == "__main__":
    # quick example-run that exercises a few functions and saves an image
    plotter = ConstructionPlotter()
    fig, ax = plotter.setup_geometry_plot()

    A = Point("A", 0.0, 0.0)
    B = Point("B", 2.0, 3.0)
    C = Point("C", 0.0, 5.0)
    D = Point("D", 1.0, 0.0)

    # plot the points and some optional connecting lines so they're visible
    plotter.plot_newpoint(ax, A)
    plotter.plot_newpoint(ax, B)
    plotter.plot_newpoint(ax, C)
    plotter.plot_newpoint(ax, D)
    prob = Problem("example", [A, B, C, D], assumptions=[], goals=[])
    dd = DDWithAR(prob)

    print("Calling find_intersection(A,B,C,D)...")
    X1 = plotter.find_intersection(A, B, C, D)
    plotter.plot_newpoint(ax, X1)
    print("Intersection:", X1, "coords:", X1.x, X1.y)

    print("Calling plot_midpoint for A,B...")
    M = plotter.plot_midpoint(ax, "M_AB", A, B)
    print("Midpoint:", M, "coords:", M.x, M.y)

    print("Drawing angle bisector for angle ABC...")
    X2 = plotter.plot_anglebisector(ax, dd, A, B, C)
    print("Angle bisector:", X2)

    I = plotter.incenter(ax, dd, A, B, C)
    O = plotter.circumcenter(ax, dd, A, B, C)

    # show figure
    plt.show()