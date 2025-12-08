import json
import random

from Read_in_Relations import *
from Read_in_Geogebra_File import *
from Constructions import *
from dd import *
from dd_ar import *
from Problem import *
from matplotlib import pyplot as plt

class DataGenerator:
    def __init__(self, points: List[RelationNode],
                 points_dict: Dict[str, Tuple[float, float]], 
                 lines: Dict[str, Tuple[float, float, float]], 
                 circles: Dict[str, Tuple[float, float, float]]):
        self.original_points = points

    def generate_rabbit(self, ddar: DDWithAR, canva: Canva):        
        choices = [
            (0.3, canva.midpoint, 2),
            (0.2, canva.anglebisector, 3),
            (0.2, canva.foot, 3),
            (0.02, canva.incenter, 3),
            (0.02, canva.incenter2, 3),
            # (0.1, canva.circumcenter, 3),
            # (0.1, canva.circumcenter2, 3),
            (0.02, canva.excenter, 3),
            (0.02, canva.excenter2, 3),
            (0.1, canva.centroid, 3),
            (0.02, canva.orthocenter, 3),
            (0.02, canva.orthocenter2, 3),
            (0.02, canva.reflect, 3),
            (0.02, canva.parallel, 3),
            (0.02, canva.tangent, 2),
            (0.02, canva.tangent2, 3)
        ]
        total_prob = sum(p for p, _, _ in choices)
        if not math.isclose(total_prob, 1.0, rel_tol=1e-9):
            raise ValueError(f"Probabilities must sum to 1.0 (got {total_prob})")

        rand_val = random.random()
        cumulative = 0.0
        for p, func, num_points in choices:
            cumulative += p
            if rand_val < cumulative:
                points_used = random.sample(ddar.problem.points, num_points)
                new_points, new_point_relations = func(*points_used)

                print(f"Constructed new point: {new_points} using points {[str(pt) for pt in points_used]} and function {func.__name__}")
                print(f"New relations: {[str(rel) for rel in new_point_relations]}")

                # invalid construction i skipped
                if (new_points is None):
                    continue

                # add everything to the solver
                if type(new_points) is not list:
                    new_points = [new_points]
                for new_point in new_points:
                    ddar.add_constructed_point(new_point)
                for r in new_point_relations:
                    ddar.add_constructed_relation(r)
                break

    def deduce_everything(self, ddar: DDWithAR):
        new_relations = []

        for rule in ddar.rules:
            new_relations.extend(rule())

        for new_rel in new_relations:
            ddar.problem.add_relation(new_rel)
    
    # TODO: get this logic working
    # returns the rabbits given a relation, string of the construction step, and the roots of the relation
    def traceback(self, relation: RelationNode, solver: DDWithAR):
        rabbits_with_roots = []
        all_relation_points = set(relation.points)
        traceback_points = []
        roots = []

        # traceback everything starting from this given relation
        def dfs(node: RelationNode):
            for point in node.points:
                if traceback_points.count(point) == 0:
                    traceback_points.append(point)
            print(node)
            print(node.parents)
            if node.parents == None or node.parents == []:
                    roots.append(node)

            for parent in node.parents:
                dfs(parent)
                
        dfs(relation)

        # TODO: roots are not being filled for some reason
        print(f"roots", roots)
        print(f"traceback points", traceback_points)

        # for all points not in the relation's points, check if there could be a rabbit for it
        for rabbit in (solver.problem.points):
            if rabbit in self.original_points or all_relation_points:
                continue
            if rabbit in traceback_points:
                # delete all the bad nodes, then add to the main list
                root_nodes = roots.copy()
                for root in root_nodes:
                    if rabbit in root.points:
                        root_nodes.remove(root)
                if root_nodes == []:
                    continue
                rabbits_with_roots.append((rabbit, root_nodes))

        return rabbits_with_roots
    
    def extract_relations(self, relations_dict: Dict[str, List[RelationNode]]) -> List[RelationNode]:
        all_relations = []
        for rel_list in relations_dict.values():
            all_relations.extend(rel_list)
        return all_relations

    def generate_new_data(self, desired_iterations: int, problem: Problem, canva: Canva):
        json_data = []
        solver = DDWithAR(problem)

        # adding N more points through random constructions
        for _ in range(desired_iterations):
            self.generate_rabbit(solver, canva)
        
        print(solver.problem.points)
        self.deduce_everything(solver)

        # for each relation, do traceback and spot out potential rabbits
        solver_relations = self.extract_relations(solver.problem.relations)
        print(f"Total relations after deduction: {solver_relations}")
        for relation in solver_relations:
            rabbits_with_roots = self.traceback(relation, solver)
            print(f"Relation: {relation.representation}, Rabbits found: {rabbits_with_roots}")
            for rabbit, rabbits_roots in rabbits_with_roots:                
                # TODO: double check if this is right
                set_up_predicates = "; ".join([root.representation for root in rabbits_roots])
                goal_predicate = relation.representation

                data_point = {
                    "input_text": f"Set_Up: {set_up_predicates}; Goal: {goal_predicate}",
                    "output_text": rabbit.representation # this should be code, not a representation, need to parse it
                }
                json_data.append(data_point)

        with open("training_data.json", "w") as f:
            json.dump(json_data, f, indent=2)

if __name__ == "__main__":
    problem_name = "problem1"

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

    data_generator = DataGenerator(points, points_dict, lines, circles)
    data_generator.generate_new_data(3, problem, canva)
    canva.plot()