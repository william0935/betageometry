"""Solve a geometry problem with the deductive engine, optionally with LLM help.

    python solve.py problem1                                 # symbolic only
    python solve.py problem1 --plot
    python solve.py usamo_2023_p1 --gemma gemma-finetuned-geometry
    python solve.py usamo_2023_p1 --random-rabbits           # baseline proposer

<<<<<<< HEAD
A problem named `foo` is read from `geogebra_files/foo.ggb` (the diagram) and
`text_files/foo.txt` (the assumptions and goal).
"""
=======
problem_name = "problem2"
>>>>>>> 6b5bb16886069cf2d557723cbaee5eb3172e2f42

import argparse
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("problem", nargs="?", default="problem1")
    parser.add_argument("--plot", action="store_true", help="show the diagram when done")
    parser.add_argument("--gemma", metavar="ADAPTER_DIR",
                        help="propose auxiliary points with a fine-tuned Gemma adapter")
    parser.add_argument("--model-id", default=None,
                        help="base checkpoint for --gemma (default: gemma-3-1b-pt)")
    parser.add_argument("--random-rabbits", action="store_true",
                        help="propose auxiliary points at random (baseline for --gemma)")
    parser.add_argument("--max-rounds", type=int, default=8,
                        help="auxiliary points to try before giving up")
    parser.add_argument("--candidates", type=int, default=4,
                        help="proposals to consider per round")
    parser.add_argument("--iterations", type=int, default=50,
                        help="deduction passes per round")
    parser.add_argument("--dump-tables", action="store_true",
                        help="write the AR tables to ar_tables.html")
    args = parser.parse_args()

    if not args.plot:
        import matplotlib
        matplotlib.use("Agg")

<<<<<<< HEAD
    from constructions import Canva
    from dd_ar import DDWithAR
    from problem import Problem
    from read_in_geogebra_file import parse_picture
    from read_in_relations import read_in_relations
    from relations import Point
    from search import solve
=======
# if problem_name == "usamo_2023_p1":
    # specific construction for USAMO 2023 P1
    # for i, p in enumerate(points):
    #     print(f"Point{i}: {p}: ({p.x}, {p.y})")
    # p, relations = canva.foot(points[0], points[1], points[2])
    
    # p, relations = canva.apply_llm_construction("foot(A, B, C)")
    # solver.add_constructed_point(p)
    # for r in relations:
    #     solver.add_constructed_relation(r)
    #     print(f"Constructed relation: {r}")
    # solver.apply_deduction_rules(50, canva)
>>>>>>> 6b5bb16886069cf2d557723cbaee5eb3172e2f42

    points_dict, lines, circles = parse_picture(f"{args.problem}.ggb")
    points = [Point(name, x, y) for name, (x, y) in points_dict.items()]
    canva = Canva(points, points_dict, lines, circles)

    assumptions, goals = read_in_relations(f"{args.problem}.txt", points)
    problem = Problem(args.problem, points, assumptions, goals)

    proposer = None
    if args.gemma:
        from gemma import DEFAULT_MODEL_ID, GemmaProposer
        proposer = GemmaProposer(model_id=args.model_id or DEFAULT_MODEL_ID,
                                 adapter_dir=args.gemma)
        print(f"Loading Gemma from {args.gemma} ...")
        proposer.load()
    elif args.random_rabbits:
        from rabbits import RandomProposer
        proposer = RandomProposer()

    solver = DDWithAR(problem, dump_tables=args.dump_tables)

    start = time.perf_counter()
    result = solve(problem, canva, proposer=proposer, solver=solver,
                   max_rounds=args.max_rounds,
                   candidates_per_round=args.candidates,
                   deduction_iterations=args.iterations)
    elapsed = time.perf_counter() - start

    print(problem)
    if result.constructions:
        print("Auxiliary constructions used:")
        for call in result.constructions:
            print(f"  {call}")
    print(f"Solved: {result.solved}")
    print(f"Time taken: {elapsed:.4f} seconds")

    if args.plot:
        canva.plot()


if __name__ == "__main__":
    main()
