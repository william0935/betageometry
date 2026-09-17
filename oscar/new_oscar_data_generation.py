import json
import random
import os
import numpy as np
import math
import sys
from itertools import combinations
from typing import List

# --- CLUSTER FIX: Must be before pyplot import ---
import matplotlib
matplotlib.use('Agg') 
from matplotlib import pyplot as plt
# -----------------------------------------------

from Read_in_Relations import *
from Read_in_Geogebra_File import *
from Constructions import *
from dd import *
from dd_ar import *
from Problem import *

# Add a single global seed so all randomness in this module is deterministic.
RANDOM_SEED = 0
random.seed(RANDOM_SEED)

def generate_rabbit(ddar: DDWithAR, canva: Canva):        
    choices = [
        (0.15, canva.midpoint, 2),
        (0.15, canva.anglebisector, 3),
        (0.15, canva.foot, 3),
        (0.05, canva.incenter, 3),
        (0.15, canva.circle, 3),
        (0.05, canva.excenter, 3),
        (0.05, canva.centroid, 3),
        (0.05, canva.orthocenter, 3),
        (0.05, canva.reflect, 3),
        (0.05, canva.parallel, 3),
        (0.05, canva.tangent, 2),
        (0.05, canva.tangent2, 3)
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

            if (new_points is None) or new_points == []:
                return generate_rabbit(ddar, canva)
            
            if type(new_points) is not list:
                new_points = [new_points]
            
            for new_point in new_points:
                for old_point in ddar.problem.points:
                    if math.isclose(new_point.x, old_point.x, rel_tol=1e-9) and \
                       math.isclose(new_point.y, old_point.y, rel_tol=1e-9):
                        return generate_rabbit(ddar, canva)

            for new_point in new_points:
                ddar.add_constructed_point(new_point)
            for r in new_point_relations:
                ddar.add_constructed_relation(r)
            return new_points, new_point_relations, [str(pt) for pt in points_used], func.__name__
    
    return generate_rabbit(ddar, canva)

def trace_back(r: RelationNode) -> List[RelationNode]:
    if r.parents == set():
        return [r]
    
    initial_relations = []
    for p in r.parents:
        p_relations = trace_back(p)
        for p_rel in p_relations:
            is_new = True
            for r in initial_relations:
                if p_rel.index == r.index:
                    is_new = False
                    break
            if is_new:
                initial_relations.append(p_rel)
    return initial_relations

def add_all_goals(solver: DDWithAR, new_points: List[Point]):
    goals = []
    for p1, p2 in combinations(solver.problem.points, 2):
        for p3, p4 in combinations(solver.problem.points, 2):
            if set({p1, p2}) != set({p3, p4}) and \
                (p1 in new_points or p2 in new_points or p3 in new_points or p4 in new_points):
                A = np.array([p1.x, p1.y])
                B = np.array([p2.x, p2.y])
                C = np.array([p3.x, p3.y])
                D = np.array([p4.x, p4.y])
                dist1 = np.linalg.norm(A - B)
                dist2 = np.linalg.norm(C - D)
                if math.isclose(dist1, dist2, rel_tol=1e-5):
                    goals.append(Congruent(p1, p2, p3, p4))
                
                grad1 = None
                if not math.isclose(B[0], A[0], rel_tol=1e-9):
                    grad1 = (B[1] - A[1]) / (B[0] - A[0])
                grad2 = None
                if not math.isclose(D[0], C[0], rel_tol=1e-9):
                    grad2 = (D[1] - C[1]) / (D[0] - C[0])
                if grad1 is None and grad2 is None:
                    goals.append(Parallel(p1, p2, p3, p4))
                elif grad1 is not None and grad2 is not None and math.isclose(grad1, grad2, rel_tol=1e-5):
                    goals.append(Parallel(p1, p2, p3, p4))

                if grad1 is None and math.isclose(grad2, 0.0, rel_tol=1e-9):
                    goals.append(Perpendicular(p1, p2, p3, p4))
                elif grad2 is None and math.isclose(grad1, 0.0, rel_tol=1e-9):
                    goals.append(Perpendicular(p1, p2, p3, p4))
                elif grad1 is not None and grad2 is not None and math.isclose(grad1 * grad2, -1.0, rel_tol=1e-5):
                    goals.append(Perpendicular(p1, p2, p3, p4))

    for p1, p2, p3 in combinations(solver.problem.points, 3):
        if (p1 in new_points or p2 in new_points or p3 in new_points):
            A = np.array([p1.x, p1.y])
            B = np.array([p2.x, p2.y])
            C = np.array([p3.x, p3.y])
            mid = (A + B) / 2
            if math.isclose(mid[0], C[0], rel_tol=1e-5) and math.isclose(mid[1], C[1], rel_tol=1e-5):
                goals.append(Midpoint(p1, p2, p3))
        
            angle_ABC = solver.problem.angle_value(p1, p2, p3)
            for p4, p5, p6 in combinations(solver.problem.points, 3):
                if (p1, p2, p3) != (p4, p5, p6):
                    angle_DEF = solver.problem.angle_value(p4, p5, p6)
                    if math.isclose(angle_ABC, angle_DEF, rel_tol=1e-5) and \
                       not math.isclose(angle_ABC, 0.0, rel_tol=1e-5) and \
                       not math.isclose(angle_ABC, 180.0, rel_tol=1e-5):
                        goals.append(EqualAngle(p1, p2, p3, p4, p5, p6))
    
    solver.problem.goals.extend(goals)
    solver.problem.remaining_goals.extend(goals)

def generate_new_data(desired_iterations: int, problem: Problem, canva: Canva, seed: int = None):
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    json_data = []
    solver = DDWithAR(problem)
    solver.apply_deduction_rules(50, canva)
    num_relations = 0
    for r_type in solver.problem.relations:
        num_relations += len(solver.problem.relations[r_type])

    for nth_point in range(desired_iterations):
        new_points, new_point_relations, points_used, func_name = generate_rabbit(solver, canva)
        
        # Stop if too many points to avoid complexity explosion
        if len(solver.problem.points) >= 16:
            print("Reached maximum number of points (16); stopping generation.")
            break
            
        func_call = func_name + "(" + ", ".join(points_used) + ")"
        add_all_goals(solver, new_points)
        solver.apply_deduction_rules(100, canva)
        
        if num_relations == len(solver.problem.relations):
            continue
        
        new_relations = []
        for r_type in solver.problem.relations:
            rels = solver.problem.relations[r_type]
            for r in rels:
                if r.index > num_relations:
                    new_relations.append(r)

        for r in new_relations:
            rabbit_check = True
            for np_ in new_points:
                if np_ in r.points:
                    rabbit_check = False
                    break
            if not rabbit_check:
                continue

            trace_back_relations = trace_back(r)
            initial_relations = []
            for tb in trace_back_relations:
                if tb in new_point_relations:
                    continue
                initial_relations.append(tb)

            input_text = "; ".join([ir.representation for ir in initial_relations]) + f"; ? {r.representation}"
            output_text = func_call
            data_point = {
                "input_text": input_text,
                "output_text": output_text
            }
            json_data.append(data_point)

        num_relations += len(new_relations)

    # Save data
    os.makedirs("data", exist_ok=True)
    filename = f"data/training_data_{seed if seed is not None else 'noset'}.json"
    with open(filename, "w") as f:
        json.dump(json_data, f, indent=2)

    return solver, json_data


if __name__ == "__main__":
    # --- MODIFIED ARGUMENT PARSING ---
    # Now accepts a filename containing a list of seeds
    if len(sys.argv) < 2:
        print("Usage: python new_oscar_data_generation.py <seed_file_path>")
        sys.exit(1)

    seed_file_path = sys.argv[1]
    
    if not os.path.exists(seed_file_path):
        print(f"Error: Seed file '{seed_file_path}' not found.")
        sys.exit(1)

    print(f"Reading seeds from {seed_file_path}...")
    
    seeds_to_process = []
    with open(seed_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.isdigit():
                seeds_to_process.append(int(line))

    print(f"Found {len(seeds_to_process)} seeds to process.")
    
    nontrivial_seeds = []
    
    for seed in seeds_to_process:
        try:
            print(f"\n=== Running seed {seed} ===")
            # Essential: Reset random state for every seed loop
            random.seed(seed)
            np.random.seed(seed)

            # Setup Problem (Canvas) for this seed
            canva = Canva([], {}, {}, {})
            points = []
            for _ in range(4):
                points.append(canva.free())

            goals = [Collinear(points[0], points[1], points[2])] # Impossible goal to drive deduction

            problem = Problem(f"random_problem_seed_{seed}", points, [], goals)

            # Run generation
            solver, json_data = generate_new_data(4, problem, canva, seed=seed)

            if json_data:
                print(f"Seed {seed} produced {len(json_data)} data points -> non-trivial")
                nontrivial_seeds.append(seed)
            else:
                print(f"Seed {seed} produced no data points")
        except Exception as e:
            print(f"Error processing seed {seed}: {e}")

    # Save summary for this specific batch file
    if nontrivial_seeds:
        os.makedirs("data", exist_ok=True)
        # Use the batch filename to create a unique summary filename
        # e.g., seed_batches/batch_0.txt -> data/nontrivial_batch_0.txt
        batch_name = os.path.basename(seed_file_path).replace('.txt', '')
        summary_filename = f"data/nontrivial_{batch_name}.txt"
        with open(summary_filename, "w") as f:
            for s in nontrivial_seeds:
                f.write(f"{s}\n")