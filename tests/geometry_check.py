"""Check a relation against the coordinates it was derived from.

The engine reasons symbolically, so a rule with a transposed index or a construction that
computes the wrong point produces relations that are *derivable* but not *true*. Nothing
in the symbolic layer can catch that; measuring the diagram can.

Angles are compared as directed angles modulo pi, which is the convention the engine
works in (see `AngleTable.add_perpendicular` and `angles_look_equal`). Comparing
undirected angles instead reports every inscribed-angle relation over a cyclic
quadrilateral as false, because the two angles are supplementary rather than equal.
"""

import math

TOL = 1e-6


def _dist(p, q):
    return math.hypot(p.x - q.x, p.y - q.y)


def _directed_angle(a, b, c):
    """Angle from line ba to line bc, in [0, pi)."""
    return (math.atan2(c.y - b.y, c.x - b.x) - math.atan2(a.y - b.y, a.x - b.x)) % math.pi


def _angles_equal(t1, t2):
    d = abs(t1 - t2)
    return min(d, math.pi - d) <= TOL


def _close(x, y):
    return abs(x - y) <= TOL * max(1.0, abs(x), abs(y))


def _cross(a, b, c):
    return abs((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x))


def _circumcircle(p1, p2, p3):
    d = 2 * (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y))
    if abs(d) < 1e-12:
        return None
    s1, s2, s3 = (p.x * p.x + p.y * p.y for p in (p1, p2, p3))
    cx = (s1 * (p2.y - p3.y) + s2 * (p3.y - p1.y) + s3 * (p1.y - p2.y)) / d
    cy = (s1 * (p3.x - p2.x) + s2 * (p1.x - p3.x) + s3 * (p2.x - p1.x)) / d
    return cx, cy, math.hypot(cx - p1.x, cy - p1.y)


def check(relation):
    """True/False if the relation is checkable numerically, None if it is not."""
    p = relation.points
    name = relation.name

    if name == "cong":
        return _close(_dist(p[0], p[1]), _dist(p[2], p[3]))

    if name == "eqangle":
        return _angles_equal(_directed_angle(p[0], p[1], p[2]),
                             _directed_angle(p[3], p[4], p[5]))

    if name == "para":
        a = math.atan2(p[1].y - p[0].y, p[1].x - p[0].x) % math.pi
        b = math.atan2(p[3].y - p[2].y, p[3].x - p[2].x) % math.pi
        d = abs(a - b)
        return min(d, math.pi - d) <= TOL

    if name == "perp":
        a = math.atan2(p[1].y - p[0].y, p[1].x - p[0].x) % math.pi
        b = math.atan2(p[3].y - p[2].y, p[3].x - p[2].x) % math.pi
        d = abs(a - b)
        return abs(min(d, math.pi - d) - math.pi / 2) <= TOL

    if name == "col":
        scale = max(1.0, _dist(p[0], p[1]) * _dist(p[0], p[2]))
        return _cross(p[0], p[1], p[2]) <= TOL * scale

    if name == "cyclic":
        circle = _circumcircle(p[0], p[1], p[2])
        if circle is None:
            return None
        cx, cy, r = circle
        return _close(math.hypot(p[3].x - cx, p[3].y - cy), r)

    if name == "circle":
        return (_close(_dist(p[0], p[1]), _dist(p[0], p[2]))
                and _close(_dist(p[0], p[1]), _dist(p[0], p[3])))

    if name == "eqratio":
        if _dist(p[2], p[3]) < 1e-12 or _dist(p[6], p[7]) < 1e-12:
            return None
        return _close(_dist(p[0], p[1]) / _dist(p[2], p[3]),
                      _dist(p[4], p[5]) / _dist(p[6], p[7]))

    if name in ("simtri1", "simtri2"):
        a = [_dist(p[0], p[1]), _dist(p[1], p[2]), _dist(p[2], p[0])]
        b = [_dist(p[3], p[4]), _dist(p[4], p[5]), _dist(p[5], p[3])]
        if min(b) < 1e-12:
            return None
        k = a[0] / b[0]
        return all(_close(a[i] / b[i], k) for i in range(3))

    if name in ("contri1", "contri2"):
        return (_close(_dist(p[0], p[1]), _dist(p[3], p[4]))
                and _close(_dist(p[1], p[2]), _dist(p[4], p[5]))
                and _close(_dist(p[2], p[0]), _dist(p[5], p[3])))

    if name == "midp":
        return (_close(_dist(p[0], p[1]), _dist(p[0], p[2]))
                and _cross(p[0], p[1], p[2]) <= TOL * max(1.0, _dist(p[1], p[2]) ** 2))

    if name == "eqarea":
        return _close(_cross(p[0], p[1], p[2]), _cross(p[3], p[4], p[5]))

    return None


def false_relations(problem):
    """Every relation in the problem that the diagram contradicts."""
    bad = []
    for relations in problem.relations.values():
        for relation in relations:
            if check(relation) is False:
                bad.append(relation)
    return bad
