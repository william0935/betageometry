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
            self.apply_congABCD_congCDEF_congABEF,
            self.apply_paraABCD_eqangleABCDCB,
            self.apply_colABC_eqangleABDBCD,
            self.apply_cong_ABDE_congBCEF_eqangleABCDEF_sameclock_ABCDE,
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
    
    # Rule: congruence of segments
    def apply_congABCD_congCDEF_congABEF(self) -> List[RelationNode]:
        new_relations = []
        congruences = self.problem.relations.get("cong", [])
        
        for i in range(len(congruences)):
            for j in range(i + 1, len(congruences)):
                transitive_rel = self.congABCD_congCDEF_congABEF(congruences[i], congruences[j])
                if transitive_rel:
                    new_relations.append(transitive_rel)
        return new_relations

    # Rule: parallel lines -> equal angles
    def apply_paraABCD_eqangleABCDCB(self) -> List[RelationNode]:
        new_relations = []
        parallels = self.problem.relations.get("para", [])
        
        for para_rel in parallels:
            angle_rel = self.paraABCD_eqangleABCDCB(para_rel)
            if angle_rel:
                new_relations.append(angle_rel)
        return new_relations
    
    # Rule: collinear points -> equal angles
    def apply_colABC_eqangleABDBCD(self) -> List[RelationNode]:
        new_relations = []
        collinears = self.problem.relations.get("col", [])
        
        for col_rel in collinears:
            angle_rels = self.colABC_eqangleABDBCD(col_rel)
            new_relations.extend(angle_rels)
        return new_relations
    
    # Rule: Triangle congruence SAS
    def apply_cong_ABDE_congBCEF_eqangleABCDEF_sameclock_ABCDE(self) -> List[RelationNode]:
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
                            contri_rel = self.cong_ABDE_congBCEF_eqangleABCDEF_sameclock_ABCDE(cong1, cong2, eqangle, sameclock)
                            if contri_rel:
                                new_relations.append(contri_rel)
        return new_relations
    
    """
    Helper methods for specific rules starts here.
    """

    # Congruent segments
    def congABCD_congCDEF_congABEF(self, cong1: Congruent, cong2: Congruent) -> Optional[Congruent]:
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
                         parents=[cong1, cong2], rule="congABCD_congCDEF_congABEF")
    
    # Parallel lines mean equal angles
    def paraABCD_eqangleABCDCB(self, para_rel: Parallel) -> Optional[EqualAngle]:
        # Parse parallel relation
        line1, line2 = list(para_rel.relation)
        A, B = list(line1)
        C, D = list(line2)
        
        # Check that the four points are not all collinear
        # If they are collinear, then the "parallel lines" are the same line
        all_points = {A, B, C, D}
        if len(all_points) < 4:
            return None  # Duplicate points
            
        # Check if all points lie on the same line by looking for collinear relations
        collinears = self.problem.relations.get("col", [])
        for col_rel in collinears:
            col_points = set(col_rel.relation)
            if all_points.issubset(col_points):
                return None  # All four points are collinear - not truly parallel lines
        
        # Create equal angle: angle ABC = angle DCB
        return EqualAngle(A, B, C, D, C, B,
                          parents=[para_rel], rule="paraABCD_eqangleABCDCB")
    
    # Collinear points mean equal angles
    def colABC_eqangleABDBCD(self, col_rel: Collinear) -> List[RelationNode]:
        new_relations = []
        points = list(col_rel.relation)
        if len(points) == 3:
            A, B, C = points[0], points[1], points[2]

            # Adding in a point D not on line ABC for each point not on the line
            for D in self.problem.points:
                if D not in points:
                    new_relations.append(EqualAngle(A, B, D, B, C, D,
                                      parents=[col_rel], rule="colABC_eqangleABDBCD"))
        return new_relations

    def cong_ABDE_congBCEF_eqangleABCDEF_sameclock_ABCDE(self, cong1: Congruent, cong2: Congruent, 
                          eqangle: EqualAngle, sameclock: Sameclock) -> Optional[CongruentTriangle1]:
        # Extract segments from congruences
        seg1_1, seg1_2 = list(cong1.relation)  # AB ≅ DE
        seg2_1, seg2_2 = list(cong2.relation)  # BC ≅ EF
        
        # Extract points
        points1_1 = list(seg1_1)  # [A, B]
        points1_2 = list(seg1_2)  # [D, E]
        points2_1 = list(seg2_1)  # [B, C]
        points2_2 = list(seg2_2)  # [E, F]
                
        # Triangle 1: A, B, C
        triangle1_points = list(set(points1_1 + points2_1))
        # Triangle 2: D, E, F
        triangle2_points = list(set(points1_2 + points2_2))
        
        if len(triangle1_points) == 3 and len(triangle2_points) == 3:
            return CongruentTriangle1(*triangle1_points, *triangle2_points,
                                    parents=[cong1, cong2, eqangle, sameclock],
                                    rule="cong_ABDE_congBCEF_eqangleABCDEF_sameclock_ABCDE")
        return None
