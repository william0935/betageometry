import json
import random

from relations import *
from Constructions import *
from dd_ar import *

class DataGenerator:
    def __init__(self, points: List[RelationNode],
                 points_dict: Dict[str, Tuple[float, float]], 
                 lines: Dict[str, Tuple[float, float, float]], 
                 circles: Dict[str, Tuple[float, float, float]]):
        self.plotter = Canva(points, points_dict, lines, circles)

    def construct_problem(self, ddar: DDWithAR):
        # initialize all the names of the angles and segments into the angle table and ratio table
        for point1, point2 in itertools.combinations(ddar.problem.points, 2):
            segment = frozenset({point1, point2})
            if segment not in ddar.angle_table.col_id:
                ddar.angle_table.add_col(segment)
            if segment not in ddar.ratio_table.col_id:
                ddar.ratio_table.add_col(segment)
        
        # add initial rules from the problem assumptions
        for cong in ddar.problem.relations.get("cong", []):
            ddar.ratio_table.add_cong(cong)
        for eqratio in ddar.problem.relations.get("eqratio", []):
            ddar.ratio_table.add_eqratio(eqratio)
        for eqangle in ddar.problem.relations.get("eqangle", []):
            ddar.angle_table.add_eqangle(eqangle)
        for para in ddar.problem.relations.get("para", []):
            ddar.angle_table.add_parallel(para)
        for col in ddar.problem.relations.get("col", []):
            ddar.angle_table.add_collinear(col)
        for perp in ddar.problem.relations.get("perp", []):
            ddar.angle_table.add_perpendicular(perp)
        for area in ddar.problem.relations.get("eqarea", []):
            ddar.area_table.add_eqarea(area)

    # TODO: get this logic working
    def traceback(self, rabbit: Point, relation: RelationNode):
        traceback_steps = []

        if rabbit in traceback_steps somehow and rabbit not in relation:
            return True

        return False
    
    def generate_rabbit(self, problem: Problem, canva: Canva) -> Point:
        choices = [
            (0.2, canva.midpoint),
            (0.4, canva.anglebisector),
            (0.1, canva.foot),
            (0.1, canva.incenter),
            (0.1, canva.incenter2),
            (0.1, canva.circumcenter),
            (0.1, canva.circumcenter2),
            (0.1, canva.excenter),
            (0.1, canva.excenter2),
            (0.1, canva.centroid),
            (0.1, canva.orthocenter),
        ]
        total_prob = sum(p for p, _ in choices)
        if not math.isclose(total_prob, 1.0, rel_tol=1e-9):
            raise ValueError(f"Probabilities must sum to 1.0 (got {total_prob})")

        rand_val = random.random()
        cumulative = 0.0
        for p, func in choices:
            cumulative += p
            if rand_val < cumulative:
                return func(problem)


    def deduce_everything(self, max_iterations: int, canva: Canva, ddar: DDWithAR):
        all_relations = []
        for iteration in range(max_iterations):
            new_point = self.generate_rabbit()
            new_relations = []

            for rule in self.rules:
                new_relations.extend(rule())

            for new_rel in new_relations:
                ddar.problem.add_relation(new_rel)
                
        return all_relations

    def generate_new_data(self, desired_iterations: int, problem: Problem, canva: Canva):
        json_data = []
        solver = DDWithAR(problem)
        self.construct_problem(solver)

        # adding a rabbit each time
        for _ in range(desired_iterations):
            rabbit = self.generate_rabbit(problem, canva)

            # call ddar
            self.deduce_everything(50, canva, solver)

            # iterate through all the derived relations and collect data
            for relation in solver.problem.relations:
                # do traceback, look at if relation uses rabbit but doesn't need it in final predicate
                if (self.traceback(rabbit, relation)):
                    # TODO: fix this formatting...
                    set_up_predicates = " ; ".join(str(r) for r in solver.problem.assumptions)
                    goal_predicate = str(relation) # is this right??
                    data_point = {
                        "input_text": f"Set_Up: {set_up_predicates}; Goal: {goal_predicate}",
                        "output_text": construction_step
                    }
                    json_data.append(data_point)

        with open("training_data.json", "w") as f:
            json.dump(json_data, f, indent=2)