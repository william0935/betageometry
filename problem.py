## Create the Problem data structure
import math
from collections import deque
from relations import *
from typing import List, Optional, Tuple
from itertools import combinations, permutations
import numpy as np

RELATION_TYPES = ["cong", "eqangle", "para", "perp", "col", "cyclic", "eqratio",
                  "simtri1", "simtri2", "contri1", "contri2", "midp", "circle", "eqarea",
                  # Inequalities. Only reasoned about when the statement contains one --
                  # see DDWithAR.uses_inequalities.
                  "gtseg", "gteseg", "gtangle", "gteangle"]

class Problem:
    def __init__(self, name: str, points: List[Point],
                 assumptions: Optional[List[RelationNode]] = None,
                 goals: Optional[List[RelationNode]] = None):
        self.name = name
        self.points = points
        self.assumptions = assumptions if assumptions is not None else []
        self.goals = goals if goals is not None else []
        self.remaining_goals = goals.copy() if goals is not None else []
        # (name, relation) -> goals waiting on it, so discharging is O(1) rather than a
        # scan of every outstanding goal on every relation added.
        self.goal_index: dict = {}
        for g in self.remaining_goals:
            self.goal_index.setdefault((g.name, g.relation), []).append(g)
        self.solved = False
        self.relations = {r : [] for r in RELATION_TYPES}
        # Mirrors `relations` as a set of hashable keys so membership tests are O(1)
        # instead of a linear scan over every relation of that type.
        self.relation_keys = {r : set() for r in RELATION_TYPES}
        self.index_counter = 1
        self.deduction_steps = []
        self.rebuild_diagram_caches()
        for r in self.assumptions:
            r.add_index(self.index_counter)
            self.index_counter += 1

        for r in self.assumptions:
            self.add_relation(r)


    def __repr__(self):
        points_str = "Points: " + ", ".join(repr(p) for p in sorted(self.points, key=lambda p: p.name))
        assumptions_str = "Assumptions:\n" + "".join(f"{r}\n" for r in self.assumptions)
        relations_str = "Known Relations:\n"
        relation_list = []
        for r_type in RELATION_TYPES:
            relation_list.extend(self.relations[r_type])
        relations_str += "".join(f"{r}\n" for r in sorted(relation_list, key=lambda r: r.index))
        goals_str = "Goals (solved):\n" if self.solved else "Goals (unsolved):\n"
        goals_str += "".join(f"{g.representation}\n" for g in self.goals)
        procedure_str = ""
        if self.solved:
            procedure_str = self.trace_back()

        return f"Problem {self.name}:\n{points_str}\n{assumptions_str}{relations_str}{goals_str}{procedure_str}"
    
    def is_solved(self):
        return self.solved

    def rebuild_diagram_caches(self):
        """Recompute the candidate sets read off the diagram.

        Called on construction and again whenever a point is added to the problem.
        """
        self.similar_triangle_pairs = self.find_similar_triangle_pairs()
        self.cyclic_quads = self.find_cyclic_quads()
        self.collinear_triples = self.find_collinear_triples()

    def add_goal(self, goal: RelationNode):
        """Register another goal, keeping `goal_index` in step with `remaining_goals`."""
        self.goals.append(goal)
        self.remaining_goals.append(goal)
        self.goal_index.setdefault((goal.name, goal.relation), []).append(goal)
        self.solved = False

    def discharge_goals(self, relation: RelationNode) -> bool:
        """Mark any outstanding goals that `relation` settles as proved.

        This is the single owner of `remaining_goals`. Callers must not remove entries
        themselves: the solver used to do its own removal as well, which raised
        ValueError as soon as two goals shared a relation, because whichever one
        `add_relation` discharged was gone by the time the solver reached it.

        All matching goals are discharged, not just the first -- data generation proposes
        goals in bulk and duplicates are routine.
        """
        key = (relation.name, relation.relation)
        matched = self.goal_index.pop(key, None)
        if not matched:
            return False

        pending = {id(g) for g in matched}
        self.remaining_goals = [g for g in self.remaining_goals if id(g) not in pending]
        self.deduction_steps.append(relation)
        if not self.remaining_goals:
            self.solved = True
        return True

    def add_relation(self, relation: RelationNode) -> Optional[List[RelationNode]]:
        new_relations = []
        if not relation.name in RELATION_TYPES:
            raise ValueError(f"Relation type {relation.name} not recognized.")

        if self.is_relation_trivial(relation):
            return None

        if relation.relation in self.relation_keys[relation.name]:
            # Already derived, so nothing new to store -- but a goal stated after the
            # fact still needs discharging.
            self.discharge_goals(relation)
            return None

        if relation.index is None:
            relation.add_index(self.index_counter)
            self.index_counter += 1

        self.relations[relation.name].append(relation)
        self.relation_keys[relation.name].add(relation.relation)
        new_relations.append(relation)

        if relation.equivalent:
            for eq in relation.equivalent:
                new_rel = self.add_relation(eq)
                if new_rel:
                    new_relations.extend(new_rel)

        self.discharge_goals(relation)

        return new_relations

    def is_relation_trivial(self, relation: RelationNode) -> bool:
        if relation.name == "cong":
            p1, p2, p3, p4 = relation.points
            return set({p1, p2}) == set({p3, p4})
        elif relation.name == "eqangle":
            p1, p2, p3, p4, p5, p6 = relation.points
            return (frozenset({p1, p2}), frozenset({p2, p3})) == (frozenset({p4, p5}), frozenset({p5, p6}))
        elif relation.name == "contri1" or relation.name == "contri2":
            p1, p2, p3, p4, p5, p6 = relation.points
            return (p1, p2, p3) == (p4, p5, p6)
        elif relation.name in ("gtseg", "gteseg", "gtangle", "gteangle"):
            # Comparing a quantity with itself: ">=" is vacuous and ">" is false.
            # Neither belongs in the database.
            left, right = relation.relation
            return left == right
        elif relation.name == "simtri1" or relation.name == "simtri2":
            p1, p2, p3, p4, p5, p6 = relation.points
            return set((p1, p2, p3)) == set((p4, p5, p6)) or \
                   relation.relation in self.relation_keys["contri1"] or \
                   relation.relation in self.relation_keys["contri2"]
        return False
    
    def trace_back(self) -> str:
        if not self.solved:
            return ""
        
        # Walk the ancestor DAG breadth-first. `seen` is keyed on object identity so a
        # node that many children share is expanded once, not once per path.
        seen = {id(s): s for s in self.deduction_steps}
        queue = deque(self.deduction_steps)
        while queue:
            step = queue.popleft()
            for parent in step.parents:
                if id(parent) not in seen:
                    seen[id(parent)] = parent
                    queue.append(parent)

        proof_steps = sorted(seen.values(), key=lambda r: r.index)

        # Number the steps for display only. Mutating `index` here would renumber nodes
        # that are still live in `self.relations`, so the numbering is kept local and the
        # parent references are rendered through the same map.
        display = {id(step): n for n, step in enumerate(proof_steps, start=1)}

        steps_str = "Deduction Steps (in order of discovery):\n"
        for step in proof_steps:
            line = f"[{display[id(step)]}] {step.representation}"
            if step.parents:
                parent_nums = sorted(display[id(p)] for p in step.parents if id(p) in display)
                if parent_nums:
                    line += f" from [{', '.join(str(n) for n in parent_nums)}]"
            if step.rule:
                line += f" via {step.rule}"
            steps_str += line + "\n"

        return steps_str

    def find_similar_triangle_pairs(self, tol=1e-5) -> List[Tuple[Point, Point, Point, Point, Point, Point]]:
        """find all possible similar triangle pairs among the points according to the diagram

        Two triangles are similar exactly when their angle multisets agree, so rather than
        testing every (combination, permutation) pair -- which is C(n,3)*P(n,3), i.e. n^6 --
        each triple's angles are computed once and the triples are sorted by their smallest
        angle. Only triples whose smallest angles agree to within `tol` can possibly match,
        so the candidate search is a sliding window over that sorted list.
        """
        points = self.points

        # One pass over the triples: O(n^3) angle computations instead of O(n^6).
        triples = []
        for p1, p2, p3 in combinations(points, 3):
            # angles at p1, p2 and p3 respectively
            a1 = self.angle_value(p2, p1, p3)
            a2 = self.angle_value(p1, p2, p3)
            a3 = self.angle_value(p2, p3, p1)
            if math.isnan(a1) or math.isnan(a2) or math.isnan(a3):
                continue
            triples.append(((p1, p2, p3), (a1, a2, a3), min(a1, a2, a3)))

        triples.sort(key=lambda t: t[2])

        pairs = []
        for i, (tri1, ang1, _) in enumerate(triples):
            p1, p2, p3 = tri1
            # Degeneracy guard, hoisted out of the inner loop: it only depends on the
            # first triangle. (A collinear triple has a 0 or 180 degree angle.)
            if ang1[1] < tol or 180 - ang1[1] < tol:
                continue

            # Candidates are the window of triples whose smallest angle is within tol.
            lo = i
            while lo > 0 and triples[i][2] - triples[lo - 1][2] <= tol:
                lo -= 1
            hi = i
            while hi + 1 < len(triples) and triples[hi + 1][2] - triples[i][2] <= tol:
                hi += 1

            for j in range(lo, hi + 1):
                tri2, ang2, _ = triples[j]
                # Every ordering of the second triangle, matching the original's
                # `permutations` over the second triple.
                for perm in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
                    q1, q2, q3 = tri2[perm[0]], tri2[perm[1]], tri2[perm[2]]
                    if p1 == q1 and p2 == q2 and p3 == q3:
                        continue
                    if abs(ang1[1] - ang2[perm[1]]) < tol and \
                       abs(ang1[2] - ang2[perm[2]]) < tol and \
                       abs(ang1[0] - ang2[perm[0]]) < tol:
                        pairs.append((p1, p2, p3, q1, q2, q3))
        return pairs

    def find_cyclic_quads(self, tol=1e-5) -> List[Tuple[Point, Point, Point, Point]]:
        """find all possible cyclic quadrilaterals among the points according to the diagram

        For a fixed chord (p1, p2) the inscribed angle at every other point is computed
        once, so the O(n^4) scan over (p3, p4) becomes pure arithmetic on cached values
        rather than O(n^4) trigonometry.
        """
        points = self.points
        quads = []
        for p1, p2 in combinations(points, 2):
            # inscribed[X] = angle p1-X-p2 ; vertex[X] = angle p1-p2-X (degeneracy guard)
            inscribed = {}
            vertex = {}
            for x in points:
                if x is p1 or x is p2:
                    continue
                inscribed[id(x)] = self.angle_value(p1, x, p2)
                vertex[id(x)] = self.angle_value(p1, p2, x)

            for p3, p4 in permutations(points, 2):
                if len({p1, p2, p3, p4}) < 4:
                    continue
                v = vertex[id(p3)]
                if v < tol or 180 - v < tol:
                    continue
                angle1 = inscribed[id(p3)]
                angle2 = inscribed[id(p4)]
                if abs(angle1 + angle2 - 180) < tol or abs(angle1 - angle2) < tol:
                    quads.append((p1, p2, p3, p4))
        return quads
    
    def find_collinear_triples(self, tol=1e-6) -> List[Tuple[Point, Point, Point]]:
        """find all possible collinear triples among the points according to the diagram

        Uses the same scale-invariant predicate as `DDWithAR.are_points_collinear`, so
        that a triple this returns is exactly a triple that predicate will accept. The
        rules driven off this list would otherwise silently miss candidates whenever the
        two tolerances disagreed.
        """
        points = self.points
        triples = []
        for p1, p2, p3 in combinations(points, 3):
            if points_look_collinear(p1, p2, p3, tol):
                triples.append((p1, p2, p3))
        return triples

    def angle_value(self, a: Point, b: Point, c: Point) -> float:
        """compute the angle value of angle ABC in degrees, in [0, 180]

        Uses atan2 on plain floats rather than building numpy arrays: this is the
        innermost operation of the diagram scans above, and the numpy version costs
        ~13us/call against ~0.4us here for an identical result.
        """
        bax, bay = a.x - b.x, a.y - b.y
        bcx, bcy = c.x - b.x, c.y - b.y
        if (bax == 0.0 and bay == 0.0) or (bcx == 0.0 and bcy == 0.0):
            # degenerate: matches the nan the previous 0/0 formulation produced
            return float("nan")
        diff = abs(math.atan2(bay, bax) - math.atan2(bcy, bcx))
        if diff > math.pi:
            diff = 2 * math.pi - diff
        return math.degrees(diff)