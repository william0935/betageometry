## Read in and parse the text of the query and convert to Problem object
## Problem statement: statement.txt

from relations import *

def parse_text(statement_path):
    f = open(statement_path, 'r')
    statement_str = f.read()
    f.close()
    pieces = statement_str.split(';')
    assumptions_strs = []
    goals_strs = []
    for str in pieces:
        str.strip()
        str.replace(';', '')
        if str[0] == '?':
            goals_strs.append(str[1:])
        else:
            assumptions_strs.append(str)
    return assumptions_strs, goals_strs

def get_relation(rel_str, all_points):
    rel_list = rel_str.split()
    word = rel_list[0]
    points_str = rel_list[1:]
    points = []
    for p_name in points_str:
        point = next((p for p in all_points if p.name == p_name), None)
        if point is None:
            raise ValueError(f"Point {p_name} not found in problem points.")
        
        points.append(point)

    if word == 'cong':
        relation = Congruent(points[0], points[1], points[2], points[3])
    elif word == 'eqangle':
        relation = EqualAngle(points[0], points[1], points[2], points[3], points[4], points[5])
    elif word == 'para':
        relation = Parallel(points[0], points[1], points[2], points[3])
    elif word == 'col':
        relation = Collinear(points[0], points[1], points[2])
    elif word == 'midp':
        relation = Midpoint(points[0], points[1], points[2])
    elif word == 'perp':
        relation = Perpendicular(points[0], points[1], points[2], points[3])
    elif word == 'sameclock':
        relation = Sameclock(points[0], points[1], points[2], points[3], points[4], points[5])
    elif word == 'cyclic':
        relation = Cyclic(points[0], points[1], points[2], points[3])
    elif word == 'simtri1':
        relation = SimilarTriangle1(points[0], points[1], points[2], points[3], points[4], points[5])
    elif word == 'simtri2':
        relation = SimilarTriangle2(points[0], points[1], points[2], points[3], points[4], points[5])
    elif word == 'contri1':
        relation = CongruentTriangle1(points[0], points[1], points[2], points[3], points[4], points[5])
    elif word == 'contri2':
        relation = CongruentTriangle2(points[0], points[1], points[2], points[3], points[4], points[5])
    else:
        raise ValueError(f"Relation type {word} not recognized.")
    
    return relation

def create_problem(assumptions_strs, goals_strs, points):
    assumptions = []
    goals = []
    for assumption in assumptions_strs:
        assumptions.append(get_relation(assumption, points))
    for goal in goals_strs:
        goals.append(get_relation(goal, points))

    return assumptions, goals

def read_in_relations(statement_path, points):
    assumptions_strs, goals_strs = parse_text(statement_path)
    assumptions, goals = create_problem(assumptions_strs, goals_strs, points)
    return assumptions, goals


def test():
    a = Point("A", 1.0, 2.0)
    b = Point("B", -1.0, -1.0)
    c = Point("C", 2.0, -2.0)
    d = Point("D", 0.5, -1.5)
    e = Point("E", 1.5, 0)
    f = Point("F", 0, 0.5)
    g = Point("G", (a.x + b.x + c.x)/3, (a.y + b.y + c.y)/3)
    points = [a, b, c, d, e, f, g]
    assumptions, goals = read_in_relations("test_relations.txt", points)
    print("Assumptions:")
    for i, r in enumerate(assumptions):
        print(f"[{i + 1}] {r}")
    print("Goals:")
    for i, r in enumerate(goals):
        print(f"[{i + 1}] {r}")

if __name__ == "__main__":
    test()