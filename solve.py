## Solve problems using deductive database

from Read_in_Relations import *
from Read_in_Geogebra_File import *
from Constructions import *
from dd import *
from dd_ar import *
from Problem import *
from matplotlib import pyplot as plt

problem_name = "problem1"

# parse info from .ggb through Read_in_Geogebra_File.py
points_dict, lines, circles = parse_picture(f"{problem_name}.ggb")

# draw initial setup
fig, ax = setup_geometry_plot()
plot_points(ax, points_dict)
plot_lines_from_eq(ax, lines)
plot_circles_from_eq(ax, circles)


# create starting points
points = []
for point in points_dict:
    point_obj = Point(point, points_dict[point][0], points_dict[point][1])
    points.append(point_obj)

# parse assumptions from .txt through Read_in_Relations.py
assumptions, goals = read_in_relations(f"{problem_name}.txt", points)

# create the problem
problem = Problem(problem_name, points, assumptions, goals)

# call dd.py or dd_ar.py until solved
solver = DDWithAR(problem)
# solver = DeductiveDatabase(problem)
solver.apply_deduction_rules(50)
print(solver.problem)

plt.show()