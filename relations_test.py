from relations import *
from problem import *

def construct_problem():
    # Define Points
    a = Point("A", 1.0, 2.0)
    b = Point("B", -1.0, -1.0)
    c = Point("C", 2.0, -2.0)
    d = Point("D", 0.5, -1.5)
    e = Point("E", 1.5, 0)
    f = Point("F", 0, 0.5)
    g = Point("G", (a.x + b.x + c.x)/3, (a.y + b.y + c.y)/3)

    # Define Relations
    midp_fab = Midpoint(f, a, b)
    midp_ebc = Midpoint(e, b, c)
    col_agd = Collinear(a, g, d)
    col_bdc = Collinear(b, d, c)
    col_bge = Collinear(b, g, e)
    col_cgf = Collinear(c, g, f)

    goal = Midpoint(g, e, f)

    return Problem(
        name="definition of the centroid",
        points={a, b, c, d, e, f, g},
        assumptions=[midp_fab, midp_ebc, col_agd, col_bdc, col_bge, col_cgf],
        goals=[goal]
    )


def main():
    test_problem = construct_problem()
    print(test_problem)


if __name__ == "__main__":
    main()