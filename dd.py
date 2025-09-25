"""
1. Look at pairs of RelationNodes 
2. Apply deduction rules to pairs to derive new RelationNodes
3. Store new RelationNodes (which automatically triggers their 'equivalent' relations)

The heavy lifting of basic equivalences is handled by relations.py through the 
'equivalent' mechanism in each RelationNode class.
"""

import itertools
from relations import *
from Problem import Problem
from typing import List, Tuple, Optional


class DeductiveDatabase:
    """
    Applies pairwise deduction rules to RelationNodes to derive new relations.
    """
    
    def __init__(self, problem: Problem):
        self.problem = problem
        
    def apply_deduction_rules(self) -> bool:
        """
        Apply all pairwise deduction rules until no new relations can be derived.
        Returns True if the problem is solved.
        """
        max_iterations = 5
        
        for iteration in range(max_iterations):
            initial_count = self._count_all_relations()
            
            # Get all current relations as a flat list
            all_relations = self._get_all_relations()
            
            # Check every pair of relations for possible deductions
            for i in range(len(all_relations)):
                for j in range(i + 1, len(all_relations)):
                    rel1, rel2 = all_relations[i], all_relations[j]
                    
                    # Try to derive new relations from this pair
                    new_relations = self._apply_pairwise_rules(rel1, rel2)
                    
                    # Add any new relations found
                    for new_rel in new_relations:
                        self.problem.add_relation(new_rel)
            
            # Check if we made progress
            if self._count_all_relations() == initial_count:
                break  # No new relations derived, stop
                
            if self.problem.is_solved():
                return True
                
        return self.problem.is_solved()
    
    def _get_all_relations(self) -> List[RelationNode]:
        """Get all relations from the problem as a flat list."""
        all_rels = []
        for rel_type in self.problem.relations:
            all_rels.extend(self.problem.relations[rel_type])
        return all_rels
    
    def _count_all_relations(self) -> int:
        """Count total number of relations."""
        return sum(len(self.problem.relations[rel_type]) for rel_type in self.problem.relations)
    
    def _apply_pairwise_rules(self, rel1: RelationNode, rel2: RelationNode) -> List[RelationNode]:
        """
        Apply deduction rules to a pair of relations to derive new ones.
        Returns list of new RelationNodes (empty if none can be derived).
        """
        new_relations = []
        
        # Rule 1: Two midpoints sharing a vertex -> medial line (parallel)
        if rel1.name == "midp" and rel2.name == "midp":
            parallel_rel = self._medial_line_rule(rel1, rel2)
            if parallel_rel:
                new_relations.append(parallel_rel)
        
        # Rule 2: Congruence transitivity
        elif rel1.name == "cong" and rel2.name == "cong":
            transitive_rel = self._congruence_transitivity(rel1, rel2)
            if transitive_rel:
                new_relations.append(transitive_rel)
        
        # Rule 3: Parallel transitivity  
        elif rel1.name == "para" and rel2.name == "para":
            transitive_rel = self._parallel_transitivity(rel1, rel2)
            if transitive_rel:
                new_relations.append(transitive_rel)
        
        # Rule 4: Collinearity extension
        elif rel1.name == "col" and rel2.name == "col":
            extended_rels = self._collinearity_extension(rel1, rel2)
            new_relations.extend(extended_rels)
        
        # Rule 5: Parallel lines -> equal angles
        elif rel1.name == "para":
            angle_rel = self._parallel_to_equal_angle(rel1)
            if angle_rel:
                new_relations.append(angle_rel)
        elif rel2.name == "para":
            angle_rel = self._parallel_to_equal_angle(rel2)
            if angle_rel:
                new_relations.append(angle_rel)

        # Rule 6: Collinear points -> equal angles
        elif rel1.name == "col":
            angle_rel = self._collinear_to_equal_angle(rel1)
            if angle_rel:
                new_relations.append(angle_rel)
        elif rel2.name == "col":
            angle_rel = self._collinear_to_equal_angle(rel2)
            if angle_rel:
                new_relations.append(angle_rel)
        
        # Check for complex multi-relation rules
        complex_rels = self._check_complex_rules(rel1, rel2)
        new_relations.extend(complex_rels)
        
        return new_relations
    
    def _medial_line_rule(self, midp1: RelationNode, midp2: RelationNode) -> Optional[RelationNode]:
        """
        Medial line theorem: If M is midpoint of AB and N is midpoint of AC, 
        then MN || BC.
        """
        # Parse midpoint representations: "midp M A B"
        parts1 = midp1.representation.split()
        parts2 = midp2.representation.split()
        
        if len(parts1) < 4 or len(parts2) < 4:
            return None
            
        mid1_name, p1_name, p2_name = parts1[1], parts1[2], parts1[3]
        mid2_name, p3_name, p4_name = parts2[1], parts2[2], parts2[3]
        
        # Check if they share exactly one vertex
        side1 = {p1_name, p2_name}
        side2 = {p3_name, p4_name}
        shared = side1 & side2
        
        if len(shared) == 1:
            # Get point objects
            mid1 = self._get_point_by_name(mid1_name)
            mid2 = self._get_point_by_name(mid2_name)
            
            # Get the two non-shared vertices
            other_vertices = (side1 | side2) - shared
            if len(other_vertices) == 2:
                vertex_names = list(other_vertices)
                v1 = self._get_point_by_name(vertex_names[0])
                v2 = self._get_point_by_name(vertex_names[1])
                
                if all([mid1, mid2, v1, v2]):
                    return Parallel(mid1, mid2, v1, v2, [midp1, midp2], "medial_line_theorem")
        
        return None
    
    def _congruence_transitivity(self, cong1: RelationNode, cong2: RelationNode) -> Optional[RelationNode]:
        """
        If AB ≅ CD and CD ≅ EF, then AB ≅ EF.
        """
        seg1, seg2 = list(cong1.relation)
        seg3, seg4 = list(cong2.relation)
        
        # Check for shared segment
        if seg2 == seg3:  # AB ≅ CD and CD ≅ EF → AB ≅ EF
            points1 = list(seg1)
            points2 = list(seg4)
            return Congruent(points1[0], points1[1], points2[0], points2[1], 
                           [cong1, cong2], "congruence_transitivity")
        elif seg1 == seg4:  # AB ≅ CD and EF ≅ AB → CD ≅ EF  
            points1 = list(seg2)
            points2 = list(seg3)
            return Congruent(points1[0], points1[1], points2[0], points2[1],
                           [cong1, cong2], "congruence_transitivity")
        
        return None
    
    def _parallel_transitivity(self, para1: RelationNode, para2: RelationNode) -> Optional[RelationNode]:
        """
        If AB || CD and CD || EF, then AB || EF.
        """
        line1, line2 = list(para1.relation)
        line3, line4 = list(para2.relation)
        
        if line2 == line3:  # AB || CD and CD || EF → AB || EF
            p1, p2 = line1
            p3, p4 = line4
            return Parallel(p1, p2, p3, p4, [para1, para2], "parallel_transitivity")
        
        return None
    
    def _collinearity_extension(self, col1: RelationNode, col2: RelationNode) -> List[RelationNode]:
        """
        If ABC are collinear and BCD are collinear, derive additional collinear relations.
        """
        points1 = set(col1.relation)
        points2 = set(col2.relation)
        
        shared = points1 & points2
        if len(shared) >= 2:  # Share at least 2 points
            all_points = points1 | points2
            if len(all_points) == 4:  # Exactly 4 distinct points total
                # All combinations of 3 points should be collinear
                new_rels = []
                for combo in itertools.combinations(all_points, 3):
                    # Skip if this combination already exists
                    combo_set = set(combo)
                    if combo_set != points1 and combo_set != points2:
                        new_rel = Collinear(*combo, [col1, col2], "collinearity_extension")
                        new_rels.append(new_rel)
                return new_rels
        
        return []
    
    def _get_point_by_name(self, name: str) -> Optional[Point]:
        """Helper to find a point by name."""
        return next((p for p in self.problem.points if p.name == name), None)
    
    def _parallel_to_equal_angle(self, para_rel: RelationNode) -> Optional[RelationNode]:
        """
        Rule 2: para A B C D -> eqangle A B C D C B
        If AB || CD, then angle ABC = angle DCB (alternate interior angles)
        """
        # Parse parallel relation
        line1, line2 = list(para_rel.relation)
        A, B = line1
        C, D = line2
        
        # Create equal angle: angle ABC = angle DCB  
        return EqualAngle(A, B, C, D, C, B, [para_rel], "parallel_equal_angles")
    
    def _collinear_to_equal_angle(self, col_rel: RelationNode) -> Optional[RelationNode]:
        """
        Rule 3: col A B C -> eqangle A B D B C D (for any point D not on line ABC)
        If A, B, C are collinear, then angle ABD = angle CBD (straight line angles)
        """
        # For simplicity, we'll need to find a suitable point D
        # This is a placeholder - in practice you'd need more context
        points = list(col_rel.relation)
        if len(points) >= 3:
            A, B, C = points[0], points[1], points[2]
            # Find a point D not on this line (simplified approach)
            for D in self.problem.points:
                if D not in points:
                    return EqualAngle(A, B, D, B, C, D, [col_rel], "collinear_equal_angles")
        return None
    
    def _check_complex_rules(self, rel1: RelationNode, rel2: RelationNode) -> List[RelationNode]:
        """
        Check for complex rules that require multiple conditions.
        Rule 1: If we have cong AB DE, cong BC EF, eqangle ABC DEF, sameclock ABC DEF
        then contri1 A B C D E F
        """
        new_relations = []
        
        # Rule 1: Triangle congruence (SAS - Side-Angle-Side)
        if rel1.name == "cong" and rel2.name == "cong":
            contri_rel = self._check_triangle_congruence_sas(rel1, rel2)
            if contri_rel:
                new_relations.append(contri_rel)
        
        return new_relations
    
    def _check_triangle_congruence_sas(self, cong1: RelationNode, cong2: RelationNode) -> Optional[RelationNode]:
        """
        Rule 1: cong A B D E; cong B C E F; eqangle A B C D E F; sameclock A B C D E F -> contri1 A B C D E F
        Check if we have two congruent sides and can find the included angle and sameclock.
        """
        # Get all current relations to check for additional conditions
        all_relations = self._get_all_relations()
        
        # Parse the two congruence relations
        seg1_1, seg1_2 = list(cong1.relation)
        seg2_1, seg2_2 = list(cong2.relation)
        
        # Try to find a pattern like AB ≅ DE and BC ≅ EF (sharing vertex B and E)
        # This is simplified - a full implementation would check all possible vertex alignments
        
        # Extract points from segments
        points1_1 = list(seg1_1)
        points1_2 = list(seg1_2)
        points2_1 = list(seg2_1) 
        points2_2 = list(seg2_2)
        
        # Look for shared vertices that could form triangles
        # This is a simplified check - would need more robust triangle detection
        if len(set(points1_1 + points2_1)) == 3 and len(set(points1_2 + points2_2)) == 3:
            # Potential triangle pattern found
            # Check for corresponding equal angle and sameclock relations
            for rel in all_relations:
                if rel.name == "eqangle":
                    # Check if this angle relation matches our triangles
                    for rel2 in all_relations:
                        if rel2.name == "sameclock":
                            # If we find matching eqangle and sameclock, create congruent triangles
                            # This is a simplified implementation
                            triangle1 = points1_1 + [p for p in points2_1 if p not in points1_1]
                            triangle2 = points1_2 + [p for p in points2_2 if p not in points1_2]
                            
                            if len(triangle1) == 3 and len(triangle2) == 3:
                                return CongruentTriangle1(*triangle1, *triangle2, 
                                                        [cong1, cong2, rel, rel2], "SAS_congruence")
        
        return None


def solve_with_deduction(problem: Problem) -> bool:
    """
    Main solver function using deductive reasoning.
    """
    print(f"Solving '{problem.name}' using deductive reasoning...")
    
    initial_relations = sum(len(problem.relations[rt]) for rt in problem.relations)
    print(f"Initial relations: {initial_relations}")
    
    # Apply deduction
    dd = DeductiveDatabase(problem)
    solved = dd.apply_deduction_rules()
    
    final_relations = sum(len(problem.relations[rt]) for rt in problem.relations)
    print(f"Final relations: {final_relations}")
    print(f"New relations derived: {final_relations - initial_relations}")
    
    if solved:
        print(f"SOLVED: {problem.name}")
    else:
        print(f"Unsolved: {problem.name}")
        print("Remaining goals:")
        for goal in problem.goals:
            print(f"  - {goal}")
    
    return solved
