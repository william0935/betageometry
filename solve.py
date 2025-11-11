## Solve problems using deductive database

from Read_in_Relations import *
from Read_in_Geogebra_File import *
from Constructions import *
from dd import *
from dd_ar import *
from Problem import *
from matplotlib import pyplot as plt
import time, datetime

problem_name = "usamo_2023_p1_rabbit"

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
end = time.perf_counter()
print(solver.problem)
print(f"Time taken: {end - start} seconds")

# TODO: call canva to create auxilliary constructions
# TODO: implement add_constructed_point, add_constructed_relation in dd_ar that updates preprocessing and AR tables

canva.plot()