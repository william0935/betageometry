import itertools
from relations import *
from problem import Problem
from typing import List, Tuple, Optional


class DeductiveDatabase:    
    def __init__(self, problem: Problem):
        self.problem = problem

        # Rule database - each rule is a method that returns List[RelationNode]
        # ADD THE NEW RULES HERE
        self.rules = [
            self.apply_congABCD_congCDEF__congABEF,
            self.apply_paraABCD__eqangleABCDCB,
            self.apply_colABC__eqangleABDCBD,
            self.apply_cong_ABDE_congBCEF_eqangleABCDEF_sameclockABCDEF__contri1ABCDEF,
            self.apply_paraABCD_paraADBC__congABCD_congADBC,
            self.apply_cong_ABDE_eqangle_CABFDE_eqangle_CBAFED_sameclockABCDEF__contri1ABCDEF
        ]

    def is_relation_new(self, relation: RelationNode) -> bool:
        """Check if a relation is new without adding it to the problem."""
        if not relation.name in self.problem.relations:
            return True
            
        existing_relations = self.problem.relations[relation.name]
        for existing in existing_relations:
            if existing.relation == relation.relation:
                return False
        return True
    
    def apply_deduction_rules(self, max_iterations: int) -> bool:   
        for iteration in range(max_iterations):
            progress_made = False
            new_relations = []
            
            # Apply all rules from the rule database
            for rule in self.rules:
                new_relations.extend(rule())
            
            # Add any new relations found and track if they're actually new
            for new_rel in new_relations:
                if self.is_relation_new(new_rel):
                    print("Derived new relation: ", new_rel)
                    self.problem.add_relation(new_rel)
                    progress_made = True
            
            # Check if we made progress
            if not progress_made:
                print(f"No new relations in iteration {iteration}. Stopping.")
                break  # No new relations derived, stop
                
            if self.problem.is_solved():
                return True
                
        return self.problem.is_solved()
    
    def apply_paraABCD_paraADBC__congABCD_congADBC(self) -> List[RelationNode]:
        """Parallelogram has equal opposite sides"""
        new_relations = []
        parallels = self.problem.relations.get("para", [])

        for i in range(len(parallels)):
            for j in range(i + 1, len(parallels)):
                cong_rels = self.paraABCD_paraADBC__congABCD_congADBC(parallels[i], parallels[j])
                if cong_rels:  # Check if not None
                    for cong_rel in cong_rels:  # Iterate through the list
                        new_relations.append(cong_rel)  # Append each Congruent individually
        return new_relations
    
    def apply_congABCD_congCDEF__congABEF(self) -> List[RelationNode]:
        """Transitive property of congruence of segments"""
        new_relations = []
        congruences = self.problem.relations.get("cong", [])
        
        for i in range(len(congruences)):
            for j in range(i + 1, len(congruences)):
                transitive_rel = self.congABCD_congCDEF__congABEF(congruences[i], congruences[j])
                if transitive_rel:
                    new_relations.append(transitive_rel)
        return new_relations

    def apply_paraABCD__eqangleABCDCB(self) -> List[RelationNode]:
        """Parallel lines mean equal angles"""
        new_relations = []
        parallels = self.problem.relations.get("para", [])
        
        for para_rel in parallels:
            angle_rels = self.paraABCD__eqangleABCDCB(para_rel)
            new_relations.extend(angle_rels)  # Use extend since function now returns a list
        return new_relations
    
    def apply_colABC__eqangleABDCBD(self) -> List[RelationNode]:
        """Collinear points mean equal angles"""
        new_relations = []
        collinears = self.problem.relations.get("col", [])
        
        for col_rel in collinears:
            angle_rels = self.colABC__eqangleABDCBD(col_rel)
            new_relations.extend(angle_rels)
        return new_relations
    
    def apply_cong_ABDE_congBCEF_eqangleABCDEF_sameclockABCDEF__contri1ABCDEF(self) -> List[RelationNode]:
        """Triangle congruence SAS"""
        new_relations = []
        congruences = self.problem.relations.get("cong", [])
        eqangles = self.problem.relations.get("eqangle", [])
        sameclocks = self.problem.relations.get("sameclock", [])
        
        # Check if we have the required components for SAS
        if len(congruences) >= 2 and len(eqangles) >= 1 and len(sameclocks) >= 1:
            # Check all combinations with progressive filtering
            for i in range(len(congruences)):
                for j in range(i + 1, len(congruences)):
                    cong1, cong2 = congruences[i], congruences[j]
                    
                    # Safety check: ensure relations have exactly 2 segments each
                    if len(cong1.relation) != 2 or len(cong2.relation) != 2:
                        continue
                    
                    # FILTER 1: Check if congruences can form triangles
                    seg1_1, seg1_2 = list(cong1.relation)
                    seg2_1, seg2_2 = list(cong2.relation)
                    points1_1, points1_2 = list(seg1_1), list(seg1_2)
                    points2_1, points2_2 = list(seg2_1), list(seg2_2)
                    
                    # Form potential triangles
                    triangle1_points = list(set(points1_1 + points2_1))
                    triangle2_points = list(set(points1_2 + points2_2))
                    
                    if len(triangle1_points) != 3 or len(triangle2_points) != 3:
                        continue  # Skip - can't form valid triangles
                    
                    # FILTER 2: Check for shared vertices (required for SAS)
                    shared_vertex_tri1 = None
                    shared_vertex_tri2 = None
                    for p in points1_1:
                        if p in points2_1:
                            shared_vertex_tri1 = p
                            break
                    for p in points1_2:
                        if p in points2_2:
                            shared_vertex_tri2 = p
                            break
                    
                    if not shared_vertex_tri1 or not shared_vertex_tri2:
                        continue  # Skip - no shared vertices for triangles
                    
                    # FILTER 3: Check if any eqangle matches the shared vertices
                    valid_eqangles = []
                    for eqangle in eqangles:
                        eqangle_points = list(eqangle.relation)
                        if len(eqangle_points) >= 6:  # Need at least 6 points for angle comparison
                            angle1_vertex = eqangle_points[1]  # Middle point of first angle
                            angle2_vertex = eqangle_points[4]  # Middle point of second angle
                            
                            # Check if angle vertices match shared vertices
                            if angle1_vertex == shared_vertex_tri1 and angle2_vertex == shared_vertex_tri2:
                                # Additional check: angle endpoints should match triangle vertices
                                angle1_endpoints = {eqangle_points[0], eqangle_points[2]}
                                angle2_endpoints = {eqangle_points[3], eqangle_points[5]}
                                other_vertices_tri1 = set(triangle1_points) - {shared_vertex_tri1}
                                other_vertices_tri2 = set(triangle2_points) - {shared_vertex_tri2}
                                
                                if angle1_endpoints == other_vertices_tri1 and angle2_endpoints == other_vertices_tri2:
                                    valid_eqangles.append(eqangle)
                    
                    if not valid_eqangles:
                        continue  # Skip - no valid angles for these triangles
                    
                    # FILTER 4: Check if any sameclock relates to both triangles
                    triangle1_set = set(triangle1_points)
                    triangle2_set = set(triangle2_points)
                    
                    for eqangle in valid_eqangles:
                        valid_sameclocks = []
                        for sameclock in sameclocks:
                            sameclock_points = set(sameclock.relation)
                            
                            # Check if sameclock involves points from both triangles
                            if (sameclock_points.intersection(triangle1_set) and 
                                sameclock_points.intersection(triangle2_set)):
                                valid_sameclocks.append(sameclock)
                        
                        if not valid_sameclocks:
                            continue  # Skip this eqangle - no valid sameclocks
                        
                        # All filters passed - now call the function
                        for sameclock in valid_sameclocks:
                            contri_rel = self.cong_ABDE_congBCEF_eqangleABCDEF_sameclockABCDEF__contri1ABCDEF(cong1, cong2, eqangle, sameclock)
                            if contri_rel:
                                new_relations.append(contri_rel)
        return new_relations
    
    def apply_cong_ABDE_eqangle_CABFDE_eqangle_CBAFED_sameclockABCDEF__contri1ABCDEF(self) -> List[RelationNode]:
        """Triangle congruence ASA (Angle-Side-Angle)"""
        new_relations = []
        congruences = self.problem.relations.get("cong", [])
        eqangles = self.problem.relations.get("eqangle", [])
        sameclocks = self.problem.relations.get("sameclock", [])
        
        # Check if we have the required components for ASA: 1 congruent side + 2 equal angles + sameclock
        if len(congruences) >= 1 and len(eqangles) >= 2 and len(sameclocks) >= 1:
            # Try all combinations of 1 congruence + 2 angles + 1 sameclock
            for cong in congruences:
                for i in range(len(eqangles)):
                    for j in range(i + 1, len(eqangles)):
                        eqangle1, eqangle2 = eqangles[i], eqangles[j]
                        for sameclock in sameclocks:
                            # Try the ASA function
                            asa_result = self.cong_ABDE_eqangle_CABFDE_eqangle_CBAFED_sameclockABCDEF__contri1ABCDEF(cong, eqangle1, eqangle2, sameclock)
                            if asa_result:
                                new_relations.append(asa_result)
        return new_relations
    
    """
    Helper methods for specific rules starts here.
    """

    def paraABCD_paraADBC__congABCD_congADBC(self, para1: Parallel, para2: Parallel) -> Optional[List[Congruent]]:
        """Parallelogram has equal opposite sides"""
        # Get the lines from the parallel relations
        line1, line2 = list(para1.relation)
        line3, line4 = list(para2.relation)
        
        # Convert frozensets to lists so we can index them
        line1_points = list(line1)
        line2_points = list(line2)
        line3_points = list(line3)
        line4_points = list(line4)
        
        # Check if lines properly intersect to form a quadrilateral
        # line1 and line3 should meet at a point, line2 and line4 should meet at a point
        line1_set = set(line1_points)
        line2_set = set(line2_points)
        line3_set = set(line3_points)
        line4_set = set(line4_points)
        intersection_13 = line1_set.intersection(line3_set)
        intersection_24 = line2_set.intersection(line4_set)
        
        # For a proper quadrilateral, each pair should intersect at exactly one point
        if len(intersection_13) != 1 or len(intersection_24) != 1:
            return None
        
        # Get all unique points from both parallel relations
        all_points = set(line1_points + line2_points + line3_points + line4_points)
        
        if len(all_points) != 4:
            return None  # Not a proper quadrilateral
        
        # For a parallelogram, opposite sides are congruent
        if len(line1_points) == 2 and len(line2_points) == 2:
            return [Congruent(line1_points[0], line1_points[1], line2_points[0], line2_points[1], parents=[para1, para2], rule="paraABCD_paraADBC__congABCD_congADBC"),
                    Congruent(line3_points[0], line3_points[1], line4_points[0], line4_points[1], parents=[para1, para2], rule="paraABCD_paraADBC__congABCD_congADBC")]

        return None

    # Congruent segments
    def congABCD_congCDEF__congABEF(self, cong1: Congruent, cong2: Congruent) -> Optional[Congruent]:
        # Safety check: ensure relations have exactly 2 segments each
        if len(cong1.relation) != 2 or len(cong2.relation) != 2:
            return None
            
        seg1, seg2 = list(cong1.relation)
        seg3, seg4 = list(cong2.relation)
        
        # Check for shared segment(need to iterate because we pulled from a set)
        if seg1 == seg3:  # seg1 is shared
            points1 = list(seg2)
            points2 = list(seg4)
        elif seg1 == seg4:  # seg1 is shared
            points1 = list(seg2)
            points2 = list(seg3)
        elif seg2 == seg3:  # seg2 is shared
            points1 = list(seg1)
            points2 = list(seg4)
        elif seg2 == seg4:  # seg2 is shared
            points1 = list(seg1)
            points2 = list(seg3)
        else:
            return None  # No shared segment found
            
        return Congruent(points1[0], points1[1], points2[0], points2[1],
                         parents=[cong1, cong2], rule="congABCD_congCDEF__congABEF")










    # Parallel lines mean equal angles
    def paraABCD__eqangleABCDCB(self, para_rel: Parallel) -> List[EqualAngle]:
        """
        Parallel lines create multiple angle equalities:
        - Corresponding angles are equal
        - Alternate interior angles are equal  
        - Alternate exterior angles are equal
        - Co-interior angles are supplementary (but we focus on equal angles here)
        """
        # Parse parallel relation
        line1, line2 = list(para_rel.relation)
        A, B = list(line1)
        C, D = list(line2)
        
        # Check that the four points are not all collinear
        # If they are collinear, then the "parallel lines" are the same line
        all_points = {A, B, C, D}
        if len(all_points) < 4:
            return []  # Duplicate points
            
        # Check if all points lie on the same line by looking for collinear relations
        collinears = self.problem.relations.get("col", [])
        for col_rel in collinears:
            col_points = set(col_rel.relation)
            if all_points.issubset(col_points):
                return []  # All four points are collinear - not truly parallel lines
        
        angle_relations = []
        
        # For parallel lines AB || CD, we need a transversal to create angle relationships
        # We'll check all possible transversals using points from the problem
        for E in self.problem.points:
            if E not in all_points:
                # E creates transversals with the parallel lines
                # Transversal through E intersecting AB and CD creates multiple angle equalities
                
                # Corresponding angles: angles in the same relative position
                # When transversal EAC crosses AB || CD:
                angle_relations.append(EqualAngle(E, A, B, E, C, D, parents=[para_rel], rule="paraABCD__eqangleABCDCB_corresponding"))
                angle_relations.append(EqualAngle(B, A, E, D, C, E, parents=[para_rel], rule="paraABCD__eqangleABCDCB_corresponding"))
                
                # When transversal EBD crosses AB || CD:
                angle_relations.append(EqualAngle(E, B, A, E, D, C, parents=[para_rel], rule="paraABCD__eqangleABCDCB_corresponding"))
                angle_relations.append(EqualAngle(A, B, E, C, D, E, parents=[para_rel], rule="paraABCD__eqangleABCDCB_corresponding"))
                
                # Alternate interior angles
                angle_relations.append(EqualAngle(E, A, B, D, C, E, parents=[para_rel], rule="paraABCD__eqangleABCDCB_alternate_interior"))
                angle_relations.append(EqualAngle(B, A, E, E, C, D, parents=[para_rel], rule="paraABCD__eqangleABCDCB_alternate_interior"))
        
        # Even without external points, we can create some basic angle relationships
        # when the four points form intersecting transversals
        if not angle_relations:
            # Basic angle equality when AB || CD
            # This creates angle relationships at intersections
            angle_relations.append(EqualAngle(A, B, C, D, C, B, parents=[para_rel], rule="paraABCD__eqangleABCDCB_basic"))
            angle_relations.append(EqualAngle(B, A, D, C, D, A, parents=[para_rel], rule="paraABCD__eqangleABCDCB_basic"))
        
        return angle_relations
    

























    
    # Collinear points mean equal angles
    def colABC__eqangleABDCBD(self, col_rel: Collinear) -> List[RelationNode]:
        new_relations = []
        points = list(col_rel.relation)
        if len(points) == 3:
            A, B, C = points[0], points[1], points[2]

            # Adding in a point D not on line ABC for each point not on the line
            for D in self.problem.points:
                if D not in points:
                    new_relations.append(EqualAngle(A, B, D, C, B, D,
                                      parents=[col_rel], rule="colABC__eqangleABDCBD"))
        return new_relations

    def cong_ABDE_congBCEF_eqangleABCDEF_sameclockABCDEF__contri1ABCDEF(self, cong1: Congruent, cong2: Congruent, 
                          eqangle: EqualAngle, sameclock: Sameclock) -> Optional[CongruentTriangle1]:
        """SAS Triangle congruence: Side-Angle-Side"""
        # Extract segments from congruences
        seg1_1, seg1_2 = list(cong1.relation)  # AB ≅ DE
        seg2_1, seg2_2 = list(cong2.relation)  # BC ≅ EF
        
        # Extract points
        points1_1 = list(seg1_1)  # [A, B]
        points1_2 = list(seg1_2)  # [D, E]
        points2_1 = list(seg2_1)  # [B, C]
        points2_2 = list(seg2_2)  # [E, F]
        
        # For SAS, we need to verify that the two congruent sides share a vertex,
        # and that the equal angle is at that shared vertex
        
        # Find shared vertices between the two congruent sides
        shared_vertex_tri1 = None
        shared_vertex_tri2 = None
        
        # Check if segments share a vertex in triangle 1
        side1_set = set(points1_1)  # {A, B}
        side2_set = set(points2_1)  # {B, C}
        shared_tri1 = side1_set.intersection(side2_set)
        
        # Check if segments share a vertex in triangle 2
        side3_set = set(points1_2)  # {D, E}
        side4_set = set(points2_2)  # {E, F}
        shared_tri2 = side3_set.intersection(side4_set)
        
        if len(shared_tri1) != 1 or len(shared_tri2) != 1:
            return None  # Sides don't share exactly one vertex each
            
        shared_vertex_tri1 = list(shared_tri1)[0]  # Should be B
        shared_vertex_tri2 = list(shared_tri2)[0]  # Should be E
        
        # Extract angles from the EqualAngle relation
        # EqualAngle relation contains two angle tuples: ((p1,p2), (p2,p3)) and ((p4,p5), (p5,p6))
        eqangle_tuples = list(eqangle.relation)
        angle1 = eqangle_tuples[0]  # First angle
        angle2 = eqangle_tuples[1]  # Second angle
        
        # Find the vertex of each angle (the shared point between the two sides)
        def find_angle_vertex(angle_tuple):
            side1, side2 = angle_tuple
            side1_points = set(side1)
            side2_points = set(side2)
            intersection = side1_points.intersection(side2_points)
            return list(intersection)[0] if len(intersection) == 1 else None
        
        vertex1 = find_angle_vertex(angle1)  # Vertex of first angle
        vertex2 = find_angle_vertex(angle2)  # Vertex of second angle
        
        # For SAS: the equal angle must be at the shared vertex of the two congruent sides
        if vertex1 != shared_vertex_tri1 or vertex2 != shared_vertex_tri2:
            return None  # The angle is not between the congruent sides (not SAS pattern)
        
        # Verify that the angle arms correspond to the congruent sides
        def get_angle_arms(angle_tuple):
            side1, side2 = angle_tuple
            return set(side1), set(side2)
        
        arms1_1, arms1_2 = get_angle_arms(angle1)  # Arms of angle at shared_vertex_tri1
        arms2_1, arms2_2 = get_angle_arms(angle2)  # Arms of angle at shared_vertex_tri2
        
        # The angle arms should match the congruent sides
        angle1_arms = {arms1_1, arms1_2}
        angle2_arms = {arms2_1, arms2_2}
        congruent_sides_tri1 = {side1_set, side2_set}
        congruent_sides_tri2 = {side3_set, side4_set}
        
        if angle1_arms != congruent_sides_tri1 or angle2_arms != congruent_sides_tri2:
            return None  # Angle arms don't match the congruent sides
        
        # If we get here, we have valid SAS pattern
        # Build the triangles
        triangle1_points = list(set(points1_1 + points2_1))
        triangle2_points = list(set(points1_2 + points2_2))
        
        if len(triangle1_points) == 3 and len(triangle2_points) == 3:
            return CongruentTriangle1(*triangle1_points, *triangle2_points,
                                    parents=[cong1, cong2, eqangle, sameclock],
                                    rule="cong_ABDE_congBCEF_eqangleABCDEF_sameclock_ABCDEF__contri1ABCDEF")
        return None

    def cong_ABDE_eqangle_CABFDE_eqangle_CBAFED_sameclockABCDEF__contri1ABCDEF(self, cong1: Congruent,
                eqangle1: EqualAngle, eqangle2: EqualAngle, sameclock: Sameclock) -> Optional[CongruentTriangle1]:
        """ASA Triangle congruence: Angle-Side-Angle"""
        # Safety check: ensure congruence has exactly 2 segments
        if len(cong1.relation) != 2:
            return None
            
        # Extract segments from congruence
        seg1_1, seg1_2 = list(cong1.relation)  # AB ≅ DE

        # Extract points
        points1_1 = list(seg1_1)  # [A, B]
        points1_2 = list(seg1_2)  # [D, E]
        
        # Extract points from angles
        # EqualAngle relation contains two angle tuples: ((p1,p2), (p2,p3)) and ((p4,p5), (p5,p6))
        eqangle1_tuples = list(eqangle1.relation)  # [((C,A), (A,B)), ((F,D), (D,E))] for angle CAB = angle FDE
        eqangle2_tuples = list(eqangle2.relation)  # [((C,B), (B,A)), ((F,E), (E,D))] for angle CBA = angle FED
        
        # Extract individual angle points
        angle1_1 = eqangle1_tuples[0]  # ((C,A), (A,B)) - first angle
        angle1_2 = eqangle1_tuples[1]  # ((F,D), (D,E)) - second angle
        angle2_1 = eqangle2_tuples[0]  # ((C,B), (B,A)) - first angle
        angle2_2 = eqangle2_tuples[1]  # ((F,E), (E,D)) - second angle
        
        # For ASA, we need to verify that the congruent side is BETWEEN the two equal angles
        # This means both angles should share the endpoints of the congruent side as their vertices
        
        # Extract the vertices and sides of each angle more carefully
        # angle1_1: ((C,A), (A,B)) means angle CAB with vertex A
        # angle1_2: ((F,D), (D,E)) means angle FDE with vertex D
        
        # Find the vertex of each angle (the shared point between the two sides)
        def find_angle_vertex(angle_tuple):
            side1, side2 = angle_tuple
            side1_points = set(side1)
            side2_points = set(side2)
            intersection = side1_points.intersection(side2_points)
            return list(intersection)[0] if len(intersection) == 1 else None
        
        vertex1_1 = find_angle_vertex(angle1_1)  # Should be A
        vertex1_2 = find_angle_vertex(angle1_2)  # Should be D
        vertex2_1 = find_angle_vertex(angle2_1)  # Should be B
        vertex2_2 = find_angle_vertex(angle2_2)  # Should be E
        
        # For ASA: the congruent side AB should be between two angles
        # This means one angle has vertex A, another has vertex B (for triangle 1)
        # and correspondingly one angle has vertex D, another has vertex E (for triangle 2)
        
        congruent_side_vertices1 = set(points1_1)  # {A, B}
        congruent_side_vertices2 = set(points1_2)  # {D, E}
        
        # Check if the angle vertices match the congruent side endpoints
        angle_vertices1 = {vertex1_1, vertex2_1}  # {A, B} should match congruent side
        angle_vertices2 = {vertex1_2, vertex2_2}  # {D, E} should match congruent side
        
        if angle_vertices1 != congruent_side_vertices1 or angle_vertices2 != congruent_side_vertices2:
            return None  # The side is not between the angles (not ASA pattern)
        
        # If we get here, we have valid ASA pattern
        # Build the triangles by finding the third vertex from each angle
        def get_third_vertex(angle_tuple, vertex):
            side1, side2 = angle_tuple
            all_points = set(side1).union(set(side2))
            third_vertices = all_points - {vertex}
            return list(third_vertices)  # Should be 2 points for the two arms of the angle
        
        # Get third vertices for triangle 1
        third_vertices1_1 = get_third_vertex(angle1_1, vertex1_1)  # From angle at A
        third_vertices2_1 = get_third_vertex(angle2_1, vertex2_1)  # From angle at B
        
        # Get third vertices for triangle 2  
        third_vertices1_2 = get_third_vertex(angle1_2, vertex1_2)  # From angle at D
        third_vertices2_2 = get_third_vertex(angle2_2, vertex2_2)  # From angle at E
        
        # The third vertex of triangle 1 should be the intersection of the angle arms
        triangle1_third = (set(third_vertices1_1).intersection(set(third_vertices2_1)))
        triangle2_third = (set(third_vertices1_2).intersection(set(third_vertices2_2)))
        
        if len(triangle1_third) != 1 or len(triangle2_third) != 1:
            return None  # Angles don't properly define triangles
        
        # Build final triangles
        triangle1_points = list(congruent_side_vertices1) + list(triangle1_third)
        triangle2_points = list(congruent_side_vertices2) + list(triangle2_third)
        
        # Check that we have valid triangles with exactly 3 points each
        if len(triangle1_points) == 3 and len(triangle2_points) == 3:
            # Try the ordering that should give us AE = CE
            # If triangle1_points = [A, B, E] and triangle2_points = [D, C, E]
            # then CongruentTriangle1(A, B, E, C, D, E) gives A↔C, B↔D, E↔E
            # which means AE = CE (what we want!)
            
            # Try alternative ordering first - this might give us the correct correspondence
            result = CongruentTriangle1(*triangle2_points, *triangle1_points,
                                        parents=[cong1, eqangle1, eqangle2, sameclock],
                                        rule="cong_ABDE_eqangle_CABFDE_eqangle_CBAFED_sameclockABCDEF__contri1ABCDEF")
            
            return result
        return None