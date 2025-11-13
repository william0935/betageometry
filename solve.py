## Solve problems using deductive database

from Read_in_Relations import *
from Read_in_Geogebra_File import *
from Constructions import *
from dd import *
from dd_ar import *
from Problem import *
from matplotlib import pyplot as plt
import time, datetime

problem_name = "Yasinsky_Geometry_Olympiad_2023_VIII_p1"

# parse info from .ggb through Read_in_Geogebra_File.py
points_dict, lines, circles = parse_picture(f"{problem_name}.ggb")

# create starting points
points = []
for point in points_dict:
    point_obj = Point(point, points_dict[point][0], points_dict[point][1])
    points.append(point_obj)

# draw initial setup
canva = Canva(points, points_dict, lines, circles)

# parse assumptions from .txt through Read_in_Relations.py
assumptions, goals = read_in_relations(f"{problem_name}.txt", points)

# create the problem
problem = Problem(problem_name, points, assumptions, goals)

# call dd.py or dd_ar.py until solved
solver = DDWithAR(problem)
# solver = DeductiveDatabase(problem)
start = time.perf_counter()
solver.apply_deduction_rules(50, canva)



# MANUAL RABBIT GENERATION

# if problem_name == "usamo_2023_p1":
#     # specific construction for USAMO 2023 P1
#     for i, p in enumerate(points):
#         print(f"Point{i}: {p}: ({p.x}, {p.y})")
#     p, relations = canva.foot(points[0], points[1], points[2])
#     solver.add_constructed_point(p)
#     for r in relations:
#         solver.add_constructed_relation(r)
#         print(f"Constructed relation: {r}")
#     solver.apply_deduction_rules(50, canva)

# if problem_name == "Yasinsky_Geometry_Olympiad_2023_VIII_p1":
#     # specific construction for Yasinsky Geometry Olympiad 2023 VIII P1
#     for i, p in enumerate(points):
#         print(f"Point{i}: {p}: ({p.x}, {p.y})")
#     p, relations = canva.mirror(points[3], points[0])
#     solver.add_constructed_point(p)
#     for r in relations:
#         solver.add_constructed_relation(r)
#         print(f"Constructed relation: {r}")
#     solver.apply_deduction_rules(50, canva)

end = time.perf_counter()
print(solver.problem)
print(f"Time taken: {end - start} seconds")

canva.plot()