## Solve problems using deductive database

from read_in_relations import *
from read_in_geogebra import *
from constructions import *
from dd import *
from problem import *
from matplotlib import pyplot as plt

problem_name = "test"

# parse info from .ggb through Read_in_Geogebra_File.py
points_dict, lines, circles = parse_picture(f"{problem_name}.ggb")

# draw initial setup with calls to Constructions.py
fig, ax = plt.subplots()
plot_point(ax, points_dict)
plot_lines_from_eq(ax, lines)
plot_circles_from_eq(ax, circles)
plt.show()

# create starting points for the problem(separate from the plotting)
points = []
for point in points_dict:
    point_obj = Point(point, points_dict[point][0], points_dict[point][1])
    points.append(point_obj)

# parse assumptions from .txt through Read_in_Relations.py
assumptions, goals = read_in_relations(f"{problem_name}.txt", points)

# create the problem
problem = Problem("problem", points, assumptions, goals)

# call dd.py until solved
solver = DeductiveDatabase(problem)
solver.apply_deduction_rules(50)
print(solver.problem)

# Manual goal checking - check if any congruent relations match our goal
goal_relation = goals[0]  # cong A E E C
print(f"\nManual goal check for: {goal_relation.representation}")
print(f"Goal relation structure: {goal_relation.relation}")

for cong_rel in solver.problem.relations["cong"]:
    if cong_rel.relation == goal_relation.relation:
        print(f"✓ GOAL ACHIEVED! Found matching relation: {cong_rel}")
        break
else:
    print("✗ Goal not found with exact relation match")
    print("Let's check all congruence relations:")
    for cong_rel in solver.problem.relations["cong"]:
        print(f"  {cong_rel} -> {cong_rel.relation}")





# # Print parents of each relation
# print("\n=== RELATION PARENTS ===")
# for relation_type, relations in solver.problem.relations.items():
#     print(f"\n{relation_type.upper()} relations:")
#     for i, relation in enumerate(relations):
#         print(f"  [{i+1}] {relation}")
#         if relation.parents:
#             print(f"      Parents: {[str(parent) for parent in relation.parents]}")
#             if relation.rule:
#                 print(f"      Rule: {relation.rule}")
#         else:
#             print(f"      Parents: None (assumption)")
#         print()