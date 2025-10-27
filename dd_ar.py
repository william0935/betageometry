import math
import itertools
from relations import *
from Problem import Problem
from typing import List, Tuple, Optional
from ar import *

class DDWithAR:    
    def __init__(self, problem: Problem):
        self.problem = problem
        self.angle_table = AngleTable([])
        self.ratio_table = RatioTable([])
        
        # Cache for triangle pair checks to avoid redundant lookups
        # Key: frozenset of two frozensets of 3 points each(for triangle pairs)
        # Value: dictionary with 'congruent', 'similar', 'checked' flags
        self.triangle_pair_cache = {}

        # Rule database - each rule is a method that returns List[RelationNode]
        # ADD THE NEW RULES HERE
        self.rules = [
            self.check_triangle_congruence_all_permutations,
            self.check_triangle_similarity_all_permutations,
            self.congMAMB_colMAB__midpointMAB,
            # self.check_triangle_congruence_SSS,
            # self.check_triangle_congruence_SAS,
            # self.check_triangle_congruence_ASA,
            # self.check_triangle_similarity_AA,
            # self.check_triangle_similarity_SAS,
            # self.check_triangle_similarity_SSS,
            self.congOAOB_congOBOC__circleOABC,
            self.circleOABC_congOCOD__cyclicABCD,
            self.eqangleABCADC__cyclicABCD,
            self.circleOABC_colOAB__perpACBC,
            self.circleOABC_perpOAAX__eqangleXABACB,
            self.circleOABC_eqangleXABACB__perpOAAX,
            self.circleOABC_midpMBC__eqangleBACBOM_eqangleCABCOM,
            self.circleOABC_colMBC_eqangleBACBOM__midpMBC,
            self.eqratioDBDCABAC_colDBC__eqangleBADDAC,
            self.eqangleBADDAC_colDBC__eqratioDBDCABAC,
            self.perpABBC_congMAMC_colMAC__congAMBM,
            self.congAPBP_congAQBQ__perpABPQ
        ]

    def print_table(self, table):
        """Print a table (AngleTable or RatioTable) in a readable format."""
        if not table.rows:
            print("  (empty)")
            print()
            return
            
        # Create header row with column indices and segment names
        header_indices = "         " + "".join(f"     [{i:1d}]" for i in range(len(table.header)))
        
        # Format segment names properly
        segment_names = []
        for col in table.header:
            if isinstance(col, frozenset):
                name = "{" + ", ".join(sorted(str(x) for x in col)) + "}"
            else:
                name = str(col)
            segment_names.append(f"{name:>8}")
        header_names = "          " + "".join(segment_names)
        
        print(header_indices)
        print(header_names)
        print("  " + "-" * (len(header_indices) - 2))  # Separator line
        
        # Print each row with proper spacing
        row_index = 0
        for relation, row_list in table.rows.items():
            for row in row_list:
                row_str = "  Row " + f"{row_index}: " + "".join(f"{val:8.2f}" if abs(val) > 1e-10 else "    0.00" for val in row)
                row_str += f"  # {relation}"
                print(row_str)
                row_index += 1
        print()  # Empty line for separation

    def apply_deduction_rules(self, max_iterations: int) -> bool:
        # initialize all the names of the angles and segments into the angle table and ratio table
        for point1, point2 in itertools.combinations(self.problem.points, 2):
            segment = frozenset({point1, point2})
            if segment not in self.angle_table.col_id:
                self.angle_table.add_col(segment)
            if segment not in self.ratio_table.col_id:
                self.ratio_table.add_col(segment)
        
        # add initial rules from the problem assumptions
        for cong in self.problem.relations.get("cong", []):
            self.ratio_table.add_cong(cong)
        for eqratio in self.problem.relations.get("eqratio", []):
            self.ratio_table.add_eqratio(eqratio)
        for eqangle in self.problem.relations.get("eqangle", []):
            self.angle_table.add_eqangle(eqangle)
        for para in self.problem.relations.get("para", []):
            self.angle_table.add_parallel(para)
        for col in self.problem.relations.get("col", []):
            self.angle_table.add_collinear(col)
        for perp in self.problem.relations.get("perp", []):
            self.angle_table.add_perpendicular(perp)

        print("Initial Angle Table:")
        self.print_table(self.angle_table)
        print("Initial Ratio Table:")
        self.print_table(self.ratio_table)

        # do iterations for dd/ar
        for iteration in range(max_iterations):
            progress_made = False
            new_relations = []
            
            # Apply all rules from the rule database
            for rule in self.rules:
                new_relations.extend(rule())
            
            # Add any new relations found and track if they're actually new
            for new_rel in new_relations:
                progress = self.problem.add_relation(new_rel)
                if progress:
                    progress_made = True
                    # print("Derived new relation: ", new_rel)
                    
                    # Update AR tables with the new relation and its equivalents
                    self.update_AR_tables_with_relation(new_rel)
                    
                    # Also update AR tables with any equivalent relations
                    if hasattr(new_rel, 'equivalent') and new_rel.equivalent:
                        for equiv_rel in new_rel.equivalent:
                            self.update_AR_tables_with_relation(equiv_rel)
            
            # Check if we made progress
            if not progress_made:
                print(f"No new relations in iteration {iteration}. Stopping.")
                print("Final Angle Table:")
                self.print_table(self.angle_table)
                print("Final Ratio Table:")
                self.print_table(self.ratio_table)
                break  # No new relations derived, stop
                
            if self.problem.is_solved():
                print(f"Problem solved in iteration {iteration}!")
                print("Final Angle Table:")
                self.print_table(self.angle_table)
                print("Final Ratio Table:")
                self.print_table(self.ratio_table)
                return True
            
            for goal in self.problem.remaining_goals:
                if self.check_relation(goal):
                    self.problem.remaining_goals.remove(goal)
                    self.problem.deduction_steps.append(goal)
                    if not self.problem.remaining_goals:
                        self.problem.solved = True
                        break
                
        return self.problem.is_solved()
    
    def congMAMB_colMAB__midpointMAB(self) -> List[RelationNode]:
        new_relations = []
        congruences = self.problem.relations.get("cong", [])
        collinears = self.problem.relations.get("col", [])
        
        for cong in congruences:
            p1, p2, p3, p4 = cong.points
            if p1 == p4:
                p3, p4 = p4, p3
            elif p2 == p3:
                p1, p2 = p2, p1
            elif p2 == p4:
                p1, p2, p3, p4 = p2, p1, p4, p3
            elif p1 == p3:
                pass
            else:
                continue

            M, A, B = p1, p2, p4
            seg_MA = frozenset({M, A})
            seg_MB = frozenset({M, B})

            # Use RatioTable (length equality), not AngleTable (direction)
            if seg_MA in self.ratio_table.col_id and seg_MB in self.ratio_table.col_id:
                cong_row = [0] * self.ratio_table.table_length()
                cong_row[self.ratio_table.col_id[seg_MA]] = 1
                cong_row[self.ratio_table.col_id[seg_MB]] = -1
                
                is_spanned, parents = self.ratio_table.is_spanned(cong_row)
                if is_spanned:
                    col = next((c for c in collinears if set(c.points) == set({M, A, B})), None)
                    if col:
                        new_relations.append(Midpoint(
                            M, A, B,
                            parents=[cong, col],
                            rule="cong_MAMB_col_MAB__midpoint_MAB"
                        ))
        
        return new_relations

    # def check_triangle_congruence_SSS(self) -> List[RelationNode]:
    #     """Check for triangle congruence using SSS criterion with caching."""
    #     new_relations = []
    #     triangles = list(itertools.combinations(self.problem.points, 3))
        
    #     for i in range(len(triangles)):
    #         for j in range(i + 1, len(triangles)):
    #             tri1 = triangles[i]
    #             tri2 = triangles[j]
                
    #             if set(tri1) == set(tri2):
    #                 continue
                
    #             cache_key = frozenset([frozenset(tri1), frozenset(tri2)])
    #             if cache_key in self.triangle_pair_cache and self.triangle_pair_cache[cache_key].get('SSS_checked'):
    #                 continue
                
    #             if cache_key not in self.triangle_pair_cache:
    #                 self.triangle_pair_cache[cache_key] = {}
    #             self.triangle_pair_cache[cache_key]['SSS_checked'] = True
                
    #             for perm in itertools.permutations(tri2):
    #                 p1, p2, p3 = tri1
    #                 p4, p5, p6 = perm
    #                 side1_1, side1_2, side1_3 = frozenset({p1, p2}), frozenset({p2, p3}), frozenset({p1, p3})
    #                 side2_1, side2_2, side2_3 = frozenset({p4, p5}), frozenset({p5, p6}), frozenset({p4, p6})
                    
    #                 if (self.are_segments_congruent(side1_1, side2_1) and
    #                     self.are_segments_congruent(side1_2, side2_2) and
    #                     self.are_segments_congruent(side1_3, side2_3)):
                        
    #                     sameclock = self.is_sameclock(p1, p2, p3, p4, p5, p6)
    #                     if sameclock:
    #                         new_relations.append(CongruentTriangle1(
    #                             p1, p2, p3, p4, p5, p6,
    #                             parents=[],
    #                             rule="check_triangle_congruence_SSS_sameclock"
    #                         ))
    #                     else:
    #                         new_relations.append(CongruentTriangle2(
    #                             p1, p2, p3, p4, p5, p6,
    #                             parents=[],
    #                             rule="check_triangle_congruence_SSS_opposite"
    #                         ))
    #                     break
        
    #     return new_relations
      
    # def check_triangle_congruence_SAS(self) -> List[RelationNode]:
    #     """Check for triangle congruence using SAS criterion with caching."""
    #     relations_to_check = self.problem.relations.get("eqangle", [])  
    #     # in the future, we will only check the new relations added in the last iteration
    #     new_relations = []
        
    #     for eqangle in relations_to_check:
    #         p1, p2, p3, p4, p5, p6 = eqangle.points

    #         seg1_1 = frozenset({p1, p2})
    #         seg1_2 = frozenset({p2, p3})
    #         seg2_1 = frozenset({p4, p5})
    #         seg2_2 = frozenset({p5, p6})

    #         # contri1 condition: p1, p2, p3, p4, p5, p6
    #         if self.are_segments_congruent(seg1_1, seg2_1) and self.are_segments_congruent(seg1_2, seg2_2):
    #             new_relations.append(CongruentTriangle1(
    #                 p1, p2, p3, p4, p5, p6,
    #                 parents=[eqangle], #in the future we will add the congruence relations as parents too
    #                 rule="check_triangle_congruence_SAS_sameclock"
    #             ))
            
    #         # contri2 condition: p1, p2, p3, p6, p5, p4
    #         if self.are_segments_congruent(seg1_1, seg2_2) and self.are_segments_congruent(seg1_2, seg2_1):
    #             new_relations.append(CongruentTriangle2(
    #                 p1, p2, p3, p6, p5, p4,
    #                 parents=[eqangle], #in the future we will add the congruence relations as parents too
    #                 rule="check_triangle_congruence_SAS_opposite"
    #             ))

    #     return new_relations
    
    # def check_triangle_congruence_ASA(self) -> List[RelationNode]:
    #     """Check for triangle congruence using ASA/AAS criterion with caching."""
    #     relations_to_check = self.problem.relations.get("eqangle", [])  
    #     new_relations = []
        
    #     for eqangle in relations_to_check:
    #         p1, p2, p3, p4, p5, p6 = eqangle.points

    #         seg1_1 = frozenset({p1, p2})
    #         seg1_2 = frozenset({p2, p3})
    #         seg1_3 = frozenset({p1, p3})
    #         seg2_1 = frozenset({p4, p5})
    #         seg2_2 = frozenset({p5, p6})
    #         seg2_3 = frozenset({p4, p6})

    #         # contri1 condition: p1, p2, p3, p4, p5, p6
    #         side_relation_found = self.are_segments_congruent(seg1_1, seg2_1) or \
    #                               self.are_segments_congruent(seg1_2, seg2_2) or \
    #                               self.are_segments_congruent(seg1_3, seg2_3)
            
    #         if side_relation_found:
    #             angle1 = (seg1_2, seg1_3)
    #             angle2 = (seg2_2, seg2_3)
    #             if self.are_angles_equal(angle1, angle2):
    #                 new_relations.append(CongruentTriangle1(
    #                     p1, p2, p3, p4, p5, p6,
    #                     parents=[eqangle], #in the future we will add the congruence relations as parents too
    #                     rule="check_triangle_congruence_ASA_sameclock"
    #                 ))
            
    #         # contri2 condition: p1, p2, p3, p6, p5, p4
    #         side_relation_found = self.are_segments_congruent(seg1_1, seg2_2) or \
    #                               self.are_segments_congruent(seg1_2, seg2_1) or \
    #                               self.are_segments_congruent(seg1_3, seg2_3)
            
    #         if side_relation_found:
    #             angle1 = (seg1_2, seg1_3)
    #             angle2 = (seg2_1, seg2_3)
    #             if self.are_angles_equal(angle1, angle2):
    #                 new_relations.append(CongruentTriangle2(
    #                     p1, p2, p3, p6, p5, p4,
    #                     parents=[eqangle], #in the future we will add the congruence relations as parents too
    #                     rule="check_triangle_congruence_ASA_opposite"
    #                 ))

    #     return new_relations
    
    # def check_triangle_similarity_AA(self) -> List[RelationNode]:
    #     """Check for triangle similarity using AA."""
    #     relations_to_check = self.problem.relations.get("eqangle", [])
    #     new_relations = []

    #     for eqangle in relations_to_check:
    #         p1, p2, p3, p4, p5, p6 = eqangle.points
    #         seg1_1 = frozenset({p1, p2})
    #         seg1_2 = frozenset({p2, p3})
    #         seg1_3 = frozenset({p1, p3})
    #         seg2_1 = frozenset({p4, p5})
    #         seg2_2 = frozenset({p5, p6})
    #         seg2_3 = frozenset({p4, p6})

    #         # simtri1 condition: p1, p2, p3, p4, p5, p6
    #         angle1 = (seg1_2, seg1_3)
    #         angle2 = (seg2_2, seg2_3)
    #         if self.are_angles_equal(angle1, angle2):
    #             new_relations.append(SimilarTriangle1(
    #                 p1, p2, p3, p4, p5, p6,
    #                 parents=[eqangle],
    #                 rule="check_triangle_similarity_AA_sameclock"
    #             ))

    #         # simtri2 condition: p1, p2, p3, p6, p5, p4
    #         angle1 = (seg1_2, seg1_3)
    #         angle2 = (seg2_1, seg2_3)
    #         if self.are_angles_equal(angle1, angle2):
    #             new_relations.append(SimilarTriangle2(
    #                 p1, p2, p3, p6, p5, p4,
    #                 parents=[eqangle],
    #                 rule="check_triangle_similarity_AA_opposite"
    #             ))

    #     return new_relations
    
    # def check_triangle_similarity_SAS(self) -> List[RelationNode]:
    #     """Check for triangle similarity using SAS."""
    #     relations_to_check = self.problem.relations.get("eqangle", [])
    #     new_relations = []

    #     for eqangle in relations_to_check:
    #         p1, p2, p3, p4, p5, p6 = eqangle.points

    #         # Sides around the equal angle and the opposite sides
    #         seg1_1 = frozenset({p1, p2})
    #         seg1_2 = frozenset({p2, p3})
    #         seg2_1 = frozenset({p4, p5})
    #         seg2_2 = frozenset({p5, p6})

    #         # Same orientation mapping: (p1,p2,p3) -> (p4,p5,p6)
    #         # Require seg1_1/seg1_2 = seg2_1/seg2_2

    #         if self.are_equal_ratio(seg1_1, seg1_2, seg2_1, seg2_2):
    #             new_relations.append(SimilarTriangle1(
    #                 p1, p2, p3, p4, p5, p6,
    #                 parents=[eqangle],
    #                 rule="check_triangle_similarity_SAS_sameclock"
    #             ))

    #         # Opposite orientation mapping: (p1,p2,p3) -> (p6,p5,p4)
    #         # Require seg1_1/seg1_2 = seg2_2/seg2_1
    #         if self.are_equal_ratio(seg1_1, seg1_2, seg2_2, seg2_1):
    #             new_relations.append(SimilarTriangle2(
    #                 p1, p2, p3, p6, p5, p4,
    #                 parents=[eqangle],
    #                 rule="check_triangle_similarity_SAS_opposite"
    #             ))

    #     return new_relations
    
    # def check_triangle_similarity_SSS(self) -> List[RelationNode]:
    #     """Check for triangle similarity using SSS."""
    #     new_relations = []
    #     # This requires checking if all three side ratios are equal
    #     # Would need more sophisticated ratio checks in AR table
    #     # Skip for now as it's complex to implement
    #     return new_relations
    
    def check_triangle_congruence_all_permutations(self) -> List[RelationNode]:
        """Check for triangle congruence for all permutations of triangle vertices."""
        new_relations = []
        triangles = list(itertools.permutations(self.problem.points, 3))
        
        for i in range(len(triangles)):
            for j in range(i + 1, len(triangles)):
                p1, p2, p3 = triangles[i]
                p4, p5, p6 = triangles[j]

                are_collinear1, _ = self.are_points_collinear(p1, p2, p3)
                are_collinear2, _ = self.are_points_collinear(p4, p5, p6)
                if are_collinear1 or are_collinear2:
                    continue

                seg1_1 = frozenset({p1, p2})
                seg1_2 = frozenset({p2, p3})
                seg1_3 = frozenset({p1, p3})
                seg2_1 = frozenset({p4, p5})
                seg2_2 = frozenset({p5, p6})
                seg2_3 = frozenset({p4, p6})

                if self.is_sameclock(p1, p2, p3, p4, p5, p6):
                    angle1_1 = (seg1_1, seg1_2)
                    angle1_2 = (seg1_1, seg1_3)
                    angle1_3 = (seg1_2, seg1_3)
                    angle2_1 = (seg2_1, seg2_2)
                    angle2_2 = (seg2_1, seg2_3)
                    angle2_3 = (seg2_2, seg2_3)

                    are_angles_equal1, parents_angles1 = self.are_angles_equal(angle1_1, angle2_1)
                    are_angles_equal2, parents_angles2 = self.are_angles_equal(angle1_2, angle2_2)
                    are_angles_equal3, parents_angles3 = self.are_angles_equal(angle1_3, angle2_3)
                    are_seg_congruent1, parents_segs1 = self.are_segments_congruent(seg1_1, seg2_1)
                    are_seg_congruent2, parents_segs2 = self.are_segments_congruent(seg1_2, seg2_2)
                    are_seg_congruent3, parents_segs3 = self.are_segments_congruent(seg1_3, seg2_3)

                    if are_angles_equal1:
                        if are_seg_congruent1 and are_seg_congruent2:
                            new_relations.append(CongruentTriangle1(
                                p1, p2, p3, p4, p5, p6,
                                parents=parents_angles1 + parents_segs1 + parents_segs2,
                                rule="check_triangle_congruence_SAS_sameclock"
                            ))
                            continue
                    if are_angles_equal2:
                        if are_seg_congruent1 and are_seg_congruent3:
                            new_relations.append(CongruentTriangle1(
                                p1, p2, p3, p4, p5, p6,
                                parents=parents_angles2 + parents_segs1 + parents_segs3,
                                rule="check_triangle_congruence_SAS_sameclock"
                            ))
                            continue
                    if are_angles_equal3:
                        if are_seg_congruent2 and are_seg_congruent3:
                            new_relations.append(CongruentTriangle1(
                                p1, p2, p3, p4, p5, p6,
                                parents=parents_angles3 + parents_segs2 + parents_segs3,
                                rule="check_triangle_congruence_SAS_sameclock"
                            ))
                            continue
                    if are_seg_congruent1 and are_seg_congruent2 and are_seg_congruent3:
                        new_relations.append(CongruentTriangle1(
                            p1, p2, p3, p4, p5, p6,
                            parents=parents_segs1 + parents_segs2 + parents_segs3,
                            rule="check_triangle_congruence_SSS_sameclock"
                        ))
                else:
                    angle1_1 = (seg1_1, seg1_2)
                    angle1_2 = (seg1_1, seg1_3)
                    angle1_3 = (seg1_2, seg1_3)
                    angle2_1 = (seg2_2, seg2_1)
                    angle2_2 = (seg2_3, seg2_1)
                    angle2_3 = (seg2_3, seg2_2)

                    are_angles_equal1, parents_angles1 = self.are_angles_equal(angle1_1, angle2_1)
                    are_angles_equal2, parents_angles2 = self.are_angles_equal(angle1_2, angle2_2)
                    are_angles_equal3, parents_angles3 = self.are_angles_equal(angle1_3, angle2_3)
                    are_seg_congruent1, parents_segs1 = self.are_segments_congruent(seg1_1, seg2_1)
                    are_seg_congruent2, parents_segs2 = self.are_segments_congruent(seg1_2, seg2_2)
                    are_seg_congruent3, parents_segs3 = self.are_segments_congruent(seg1_3, seg2_3)
                    
                    if are_angles_equal1:
                        if are_seg_congruent1 and are_seg_congruent2:
                            new_relations.append(CongruentTriangle2(
                                p1, p2, p3, p4, p5, p6,
                                parents=parents_angles1 + parents_segs1 + parents_segs2,
                                rule="check_triangle_congruence_SAS_opposite"
                            ))
                            continue
                    if are_angles_equal2:
                        if are_seg_congruent1 and are_seg_congruent3:
                            new_relations.append(CongruentTriangle2(
                                p1, p2, p3, p4, p5, p6,
                                parents=parents_angles2 + parents_segs1 + parents_segs3,
                                rule="check_triangle_congruence_SAS_opposite"
                            ))
                            continue
                    if are_angles_equal3:
                        if are_seg_congruent2 and are_seg_congruent3:
                            new_relations.append(CongruentTriangle2(
                                p1, p2, p3, p4, p5, p6,
                                parents=parents_angles3 + parents_segs2 + parents_segs3,
                                rule="check_triangle_congruence_SAS_opposite"
                            ))
                            continue
                    if are_seg_congruent1 and are_seg_congruent2 and are_seg_congruent3:
                        new_relations.append(CongruentTriangle2(
                            p1, p2, p3, p4, p5, p6,
                            parents=parents_segs1 + parents_segs2 + parents_segs3,
                            rule="check_triangle_congruence_SSS_opposite"
                        ))
        
        return new_relations

    def check_triangle_similarity_all_permutations(self) -> List[RelationNode]:
        """Check for triangle similarity for all permutations of triangle vertices."""
        new_relations = []
        triangles = list(itertools.permutations(self.problem.points, 3))

        for i in range(len(triangles)):
            for j in range(i + 1, len(triangles)):
                p1, p2, p3 = triangles[i]
                p4, p5, p6 = triangles[j]

                are_collinear1, _ = self.are_points_collinear(p1, p2, p3)
                are_collinear2, _ = self.are_points_collinear(p4, p5, p6)
                if are_collinear1 or are_collinear2:
                    continue

                seg1_1 = frozenset({p1, p2})
                seg1_2 = frozenset({p2, p3})
                seg1_3 = frozenset({p1, p3})
                seg2_1 = frozenset({p4, p5})
                seg2_2 = frozenset({p5, p6})
                seg2_3 = frozenset({p4, p6})

                if self.is_sameclock(p1, p2, p3, p4, p5, p6):
                    angle1_1 = (seg1_1, seg1_2)
                    angle1_2 = (seg1_1, seg1_3)
                    angle1_3 = (seg1_2, seg1_3)
                    angle2_1 = (seg2_1, seg2_2)
                    angle2_2 = (seg2_1, seg2_3)
                    angle2_3 = (seg2_2, seg2_3)

                    are_angles_equal1, parents_angles1 = self.are_angles_equal(angle1_1, angle2_1)
                    are_angles_equal2, parents_angles2 = self.are_angles_equal(angle1_2, angle2_2)
                    are_angles_equal3, parents_angles3 = self.are_angles_equal(angle1_3, angle2_3)
                    are_equal_ratio1, parents_ratio1 = self.are_equal_ratio(seg1_1, seg1_2, seg2_1, seg2_2)
                    are_equal_ratio2, parents_ratio2 = self.are_equal_ratio(seg1_1, seg1_3, seg2_1, seg2_3)
                    are_equal_ratio3, parents_ratio3 = self.are_equal_ratio(seg1_2, seg1_3, seg2_2, seg2_3)

                    if are_angles_equal1 and are_angles_equal2:
                        new_relations.append(SimilarTriangle1(
                            p1, p2, p3, p4, p5, p6,
                            parents=parents_angles1 + parents_angles2,
                            rule="check_triangle_similarity_AA_sameclock"
                        ))
                        continue
                    if are_angles_equal1:
                        if are_equal_ratio1:
                            new_relations.append(SimilarTriangle1(
                                p1, p2, p3, p4, p5, p6,
                                parents=parents_angles1 + parents_ratio1,
                                rule="check_triangle_similarity_SAS_sameclock"
                            ))
                            continue
                    if are_angles_equal2:
                        if are_equal_ratio2:
                            new_relations.append(SimilarTriangle1(
                                p1, p2, p3, p4, p5, p6,
                                parents=parents_angles2 + parents_ratio2,
                                rule="check_triangle_similarity_SAS_sameclock"
                            ))
                            continue
                    if are_angles_equal3:
                        if are_equal_ratio3:
                            new_relations.append(SimilarTriangle1(
                                p1, p2, p3, p4, p5, p6,
                                parents=parents_angles3 + parents_ratio3,
                                rule="check_triangle_similarity_SAS_sameclock"
                            ))
                            continue
                    if are_equal_ratio1 and are_equal_ratio2:
                        new_relations.append(SimilarTriangle1(
                            p1, p2, p3, p4, p5, p6,
                            parents=parents_ratio1 + parents_ratio2,
                            rule="check_triangle_similarity_SSS_sameclock"
                        ))
                else:
                    angle1_1 = (seg1_1, seg1_2)
                    angle1_2 = (seg1_1, seg1_3)
                    angle1_3 = (seg1_2, seg1_3)
                    angle2_1 = (seg2_2, seg2_1)
                    angle2_2 = (seg2_3, seg2_1)
                    angle2_3 = (seg2_3, seg2_2)

                    are_angles_equal1, parents_angles1 = self.are_angles_equal(angle1_1, angle2_1)
                    are_angles_equal2, parents_angles2 = self.are_angles_equal(angle1_2, angle2_2)
                    are_angles_equal3, parents_angles3 = self.are_angles_equal(angle1_3, angle2_3)
                    are_equal_ratio1, parents_ratio1 = self.are_equal_ratio(seg1_1, seg1_2, seg2_1, seg2_2)
                    are_equal_ratio2, parents_ratio2 = self.are_equal_ratio(seg1_1, seg1_3, seg2_1, seg2_3)
                    are_equal_ratio3, parents_ratio3 = self.are_equal_ratio(seg1_2, seg1_3, seg2_2, seg2_3)
    
                    if are_angles_equal1 and are_angles_equal2:
                        new_relations.append(SimilarTriangle2(
                            p1, p2, p3, p4, p5, p6,
                            parents=parents_angles1 + parents_angles2,
                            rule="check_triangle_similarity_AA_opposite"
                        ))
                        continue
                    if are_angles_equal1:
                        if are_equal_ratio1:
                            new_relations.append(SimilarTriangle2(
                                p1, p2, p3, p4, p5, p6,
                                parents=parents_angles1 + parents_ratio1,
                                rule="check_triangle_similarity_SAS_opposite"
                            ))
                            continue
                    if are_angles_equal2:
                        if are_equal_ratio2:
                            new_relations.append(SimilarTriangle2(
                                p1, p2, p3, p4, p5, p6,
                                parents=parents_angles2 + parents_ratio2,
                                rule="check_triangle_similarity_SAS_opposite"
                            ))
                            continue
                    if are_angles_equal3:
                        if are_equal_ratio3:
                            new_relations.append(SimilarTriangle2(
                                p1, p2, p3, p4, p5, p6,
                                parents=parents_angles3 + parents_ratio3,
                                rule="check_triangle_similarity_SAS_opposite"
                            ))
                            continue
                    if are_equal_ratio1 and are_equal_ratio2:
                        new_relations.append(SimilarTriangle2(
                            p1, p2, p3, p4, p5, p6,
                            parents=parents_ratio1 + parents_ratio2,
                            rule="check_triangle_similarity_SSS_opposite"
                        ))

        return new_relations

    def congOAOB_congOBOC__circleOABC(self) -> List[RelationNode]:
        """Check for circle formation from congruent segments."""
        relations_to_check = self.problem.relations.get("cong", [])
        new_relations = []
        for cong1 in relations_to_check:
            p1, p2, p3, p4 = self.form_triangle(cong1.points)
            if p1 is None:
                continue

            for p in self.problem.points:
                if p in {p1, p2, p4}:
                    continue
                seg1 = frozenset({p1, p2})
                seg2 = frozenset({p1, p4})
                seg3 = frozenset({p1, p})

                if self.are_segments_congruent(seg1, seg3) or \
                   self.are_segments_congruent(seg2, seg3):
                    new_relations.append(Circle(
                        p1, p2, p4, p,
                        parents=[cong1],
                        rule="congOAOB_congOBOC__circleOABC"
                    ))

        return new_relations
    
    def circleOABC_congOCOD__cyclicABCD(self) -> List[RelationNode]:
        """Check for cyclic quadrilateral formation from circle and congruent segments."""
        relations_to_check = self.problem.relations.get("circle", [])
        new_relations = []
        for circle in relations_to_check:
            p1, p2, p3, p4 = circle.points
            for p in self.problem.points:
                if p in {p1, p2, p3, p4}:
                    continue
                seg1 = frozenset({p1, p2})
                seg2 = frozenset({p1, p3})
                seg3 = frozenset({p1, p4})
                seg4 = frozenset({p1, p})

                if self.are_segments_congruent(seg1, seg4) or \
                   self.are_segments_congruent(seg2, seg4) or \
                   self.are_segments_congruent(seg3, seg4):
                    new_relations.append(Cyclic(
                        p2, p3, p4, p,
                        parents=[circle],
                        rule="circleOABC_congOCOD__cyclicABCD"
                    ))

        return new_relations
    
    def eqangleABCADC__cyclicABCD(self) -> List[RelationNode]:
        """Check for cyclic quadrilateral formation from equal angles."""
        relations_to_check = self.problem.relations.get("eqangle", [])
        new_relations = []
        for eqangle in relations_to_check:
            p1, p2, p3, p4, p5, p6 = eqangle.points
            if p1 == p4 and p3 == p6:
                new_relations.append(Cyclic(
                    p1, p2, p3, p5,
                    parents=[eqangle],
                    rule="eqangleABCADC__cyclicABCD"
                ))

        return new_relations
    
    def circleOABC_colOAB__perpACBC(self) -> List[RelationNode]:
        """Check for perpendicularity from circle and collinearity."""
        relations_to_check = self.problem.relations.get("circle", [])
        new_relations = []
        
        for circle in relations_to_check:
            p1, p2, p3, p4 = circle.points
            if self.are_points_collinear(p1, p2, p3):
                new_relations.append(Perpendicular(
                    p2, p4, p3, p4,
                    parents=[circle],
                    rule="circleOABC_colOAB__perpACBC"
                ))
            
            if self.are_points_collinear(p1, p2, p4):
                new_relations.append(Perpendicular(
                    p2, p3, p3, p4,
                    parents=[circle],
                    rule="circleOABC_colOAC__perpABBC"
                ))

            if self.are_points_collinear(p1, p3, p4):
                new_relations.append(Perpendicular(
                    p2, p3, p2, p4,
                    parents=[circle],
                    rule="circleOABC_colOBC__perpABAC"
                ))

        return new_relations

    def circleOABC_perpOAAX__eqangleXABACB(self) -> List[RelationNode]:
        """Check for equal angles from circle and perpendicularity."""
        relations_to_check = self.problem.relations.get("circle", [])
        new_relations = []
        
        for circle in relations_to_check:
            p1, p2, p3, p4 = circle.points
            for p in self.problem.points:
                if p in {p1, p2, p3, p4}:
                    continue

                if self.are_points_perpendicular(p1, p2, p2, p):
                    new_relations.append(EqualAngle(
                        p, p2, p3, p2, p4, p3,
                        parents=[circle],
                        rule="circleOABC_perpOAAX__eqangleXABACB"
                    ))
                if self.are_points_perpendicular(p1, p3, p3, p):
                    new_relations.append(EqualAngle(
                        p, p3, p2, p3, p4, p2,
                        parents=[circle],
                        rule="circleOABC_perpOAAX__eqangleXABACB"
                    ))
                if self.are_points_perpendicular(p1, p4, p4, p):
                    new_relations.append(EqualAngle(
                        p, p4, p2, p4, p3, p2,
                        parents=[circle],
                        rule="circleOABC_perpOAAX__eqangleXABACB"
                    ))

        return new_relations
    
    def circleOABC_eqangleXABACB__perpOAAX(self) -> List[RelationNode]:
        """Check for perpendicularity from circle and equal angles."""
        relations_to_check = self.problem.relations.get("circle", [])
        new_relations = []
        
        for circle in relations_to_check:
            p1, p2, p3, p4 = circle.points
            for p in self.problem.points:
                if p in {p1, p2, p3, p4}:
                    continue

                if self.are_angles_equal((frozenset({p, p2}), frozenset({p2, p3})),
                                         (frozenset({p, p3}), frozenset({p3, p4}))):
                    new_relations.append(Perpendicular(
                        p1, p2, p2, p,
                        parents=[circle],
                        rule="circleOABC_eqangleXABACB__perpOAAX"
                    ))
                if self.are_angles_equal((frozenset({p, p3}), frozenset({p3, p2})),
                                         (frozenset({p, p2}), frozenset({p2, p4}))):
                    new_relations.append(Perpendicular(
                        p1, p3, p3, p,
                        parents=[circle],
                        rule="circleOABC_eqangleXABACB__perpOAAX"
                    ))
                if self.are_angles_equal((frozenset({p, p4}), frozenset({p4, p2})),
                                         (frozenset({p, p2}), frozenset({p2, p3}))):
                    new_relations.append(Perpendicular(
                        p1, p4, p4, p,
                        parents=[circle],
                        rule="circleOABC_eqangleXABACB__perpOAAX"
                    ))

        return new_relations

    def circleOABC_midpMBC__eqangleBACBOM_eqangleCABCOM(self) -> List[RelationNode]:
        """Check for equal angles from circle and midpoint."""
        relations_to_check = self.problem.relations.get("circle", [])
        new_relations = []
        
        for circle in relations_to_check:
            p1, p2, p3, p4 = circle.points
            for p in self.problem.points:
                if p in {p1, p2, p3, p4}:
                    continue

                is_midpoint, parents = self.is_midpoint(p, p3, p4)
                if is_midpoint:
                    new_relations.append(EqualAngle(
                        p3, p2, p4, p3, p1, p,
                        parents=[circle],
                        rule="circleOABC_midpMBC__eqangleBACBOM"
                    )).append(EqualAngle(
                        p4, p2, p3, p4, p1, p,
                        parents=[circle],
                        rule="circleOABC_midpMBC__eqangleCABCOM"
                    ))
                
                is_midpoint, parents = self.is_midpoint(p, p2, p4)
                if is_midpoint:
                    new_relations.append(EqualAngle(
                        p2, p3, p4, p2, p1, p,
                        parents=[circle],
                        rule="circleOABC_midpMBC__eqangleBACBOM"
                    )).append(EqualAngle(
                        p4, p3, p2, p4, p1, p,
                        parents=[circle],
                        rule="circleOABC_midpMBC__eqangleCABCOM"
                    ))
                
                is_midpoint, parents = self.is_midpoint(p, p2, p3)
                if is_midpoint:
                    new_relations.append(EqualAngle(
                        p2, p4, p3, p2, p1, p,
                        parents=[circle],
                        rule="circleOABC_midpMBC__eqangleBACBOM"
                    )).append(EqualAngle(
                        p3, p4, p2, p3, p1, p,
                        parents=[circle],
                        rule="circleOABC_midpMBC__eqangleCABCOM"
                    ))

        return new_relations
    
    def circleOABC_colMBC_eqangleBACBOM__midpMBC(self) -> List[RelationNode]:
        """Check for midpoint from circle, collinearity, and equal angles."""
        relations_to_check = self.problem.relations.get("circle", [])
        new_relations = []
        
        for circle in relations_to_check:
            p1, p2, p3, p4 = circle.points
            for p in self.problem.points:
                if p in {p1, p2, p3, p4}:
                    continue

                if self.are_points_collinear(p3, p4, p) and \
                   self.are_angles_equal(
                       (frozenset({p3, p2}), frozenset({p2, p4})),
                       (frozenset({p3, p1}), frozenset({p1, p}))
                   ):
                    new_relations.append(Midpoint(
                        p, p3, p4,
                        parents=[circle],
                        rule="circleOABC_colMBC_eqangleBACBOM__midpMBC"
                    ))
                if self.are_points_collinear(p2, p4, p) and \
                   self.are_angles_equal(
                       (frozenset({p2, p3}), frozenset({p3, p4})),
                       (frozenset({p2, p1}), frozenset({p1, p}))
                   ):
                    new_relations.append(Midpoint(
                        p, p2, p4,
                        parents=[circle],
                        rule="circleOABC_colMBC_eqangleCABCOM__midpMBC"
                    ))
                if self.are_points_collinear(p2, p3, p) and \
                   self.are_angles_equal(
                       (frozenset({p2, p4}), frozenset({p4, p3})),
                       (frozenset({p2, p1}), frozenset({p1, p}))
                   ):
                    new_relations.append(Midpoint(
                        p, p2, p3,
                        parents=[circle],
                        rule="circleOABC_colMBC_eqangleCABCOM__midpMBC"
                    ))

        return new_relations

    def eqratioDBDCABAC_colDBC__eqangleBADDAC(self) -> List[RelationNode]:
        """Check for equal angles from equal ratios and collinearity."""
        new_relations = []
        for p1, p2, p3, p4 in itertools.permutations(self.problem.points, 4):
            seg1 = frozenset({p4, p2})
            seg2 = frozenset({p4, p3})
            seg3 = frozenset({p1, p2})
            seg4 = frozenset({p1, p3})

            if self.are_equal_ratio(seg1, seg2, seg3, seg4) and \
               self.are_points_collinear(p2, p3, p4) and \
               not self.are_points_collinear(p1, p2, p3):
                new_relations.append(EqualAngle(
                    p2, p1, p4, p4, p1, p3,
                    parents=[],
                    rule="eqratioDBDCABAC_colDBC__eqangleBADDAC"
                ))

        return new_relations
    
    def eqangleBADDAC_colDBC__eqratioDBDCABAC(self) -> List[RelationNode]:
        """Check for equal ratios from equal angles and collinearity."""
        relations_to_check = self.problem.relations.get("eqangle", [])
        new_relations = []
        for eqangle in relations_to_check:
            p2, p1, p4, p4b, p1b, p3 = eqangle.points
            if p4 != p4b or p1 != p1b:
                continue

            if self.are_points_collinear(p2, p3, p4) and \
               not self.are_points_collinear(p1, p2, p3):
                new_relations.append(EqualRatio(
                    p4, p2, p4, p3, p1, p2, p1, p3,
                    parents=[eqangle],
                    rule="eqangleBADDAC_colDBC__eqratioDBDCABAC"
                ))

        return new_relations

    def perpABBC_congMAMC_colMAC__congAMBM(self) -> List[RelationNode]:
        """Check for congruent segments from perpendicularity and collinearity."""
        relations_to_check = self.problem.relations.get("perp", [])
        new_relations = []
        for perp in relations_to_check:
            p2, p1, _, p3 = self.form_triangle(perp.points)
            if p2 is None:
                continue

            for p in self.problem.points:
                if p in {p1, p2, p3}:
                    continue
                if self.are_segments_congruent(frozenset({p, p1}), frozenset({p, p3})) and \
                   self.are_points_collinear(p1, p3, p):
                    new_relations.append(Congruent(
                        p1, p, p2, p,
                        parents=[perp],
                        rule="perpABBC_congMAMC_colMAC__congAMBM"
                    ))

        return new_relations

    def congAPBP_congAQBQ__perpABPQ(self) -> List[RelationNode]:
        """Check for perpendicularity from congruent segments."""
        relations_to_check = self.problem.relations.get("cong", [])
        new_relations = []
        for cong1 in relations_to_check:
            p3, p1, _, p2 = self.form_triangle(cong1.points)
            if p3 is None:
                continue

            for p in self.problem.points:
                if p in {p1, p2, p3}:
                    continue
                if self.are_segments_congruent(frozenset({p1, p}), frozenset({p2, p})):
                    new_relations.append(Perpendicular(
                        p1, p2, p3, p,
                        parents=[cong1],
                        rule="congAPBP_congAQBQ__perpABPQ"
                    ))
            
        return new_relations
    

    """
    Helper methods for specific rules starts here.
    """

    def check_relation(self, rel: RelationNode) -> bool:
        """Check if a specific relation is discovered."""
        if isinstance(rel, Congruent):
            p1, p2, p3, p4 = rel.points
            seg1 = frozenset({p1, p2})
            seg2 = frozenset({p3, p4})
            return self.are_segments_congruent(seg1, seg2)
        elif isinstance(rel, EqualRatio):
            p1, p2, p3, p4, p5, p6, p7, p8 = rel.points
            seg1 = frozenset({p1, p2})
            seg2 = frozenset({p3, p4})
            seg3 = frozenset({p5, p6})
            seg4 = frozenset({p7, p8})
            return self.are_equal_ratio(seg1, seg2, seg3, seg4)
        elif isinstance(rel, EqualAngle):
            p1, p2, p3, p4, p5, p6 = rel.points
            angle1 = (frozenset({p1, p2}), frozenset({p2, p3}))
            angle2 = (frozenset({p4, p5}), frozenset({p5, p6}))
            return self.are_angles_equal(angle1, angle2)
        elif isinstance(rel, Collinear):
            p1, p2, p3 = rel.points
            return self.are_points_collinear(p1, p2, p3)
        elif isinstance(rel, Parallel):
            p1, p2, p3, p4 = rel.points
            return self.are_points_parallel(p1, p2, p3, p4)
        elif isinstance(rel, Perpendicular):
            p1, p2, p3, p4 = rel.points
            return self.are_points_perpendicular(p1, p2, p3, p4)
        else:
            return False
    
    def update_AR_tables_with_relation(self, rel: RelationNode):
        """Update AR tables when a new relation is discovered."""
        if isinstance(rel, Congruent):
            self.ratio_table.add_cong(rel)
        elif isinstance(rel, EqualRatio):
            self.ratio_table.add_eqratio(rel)
        elif isinstance(rel, EqualAngle):
            self.angle_table.add_eqangle(rel)
        elif isinstance(rel, Collinear):
            self.angle_table.add_collinear(rel)
        elif isinstance(rel, Parallel):
            self.angle_table.add_parallel(rel)
        elif isinstance(rel, Perpendicular):
            self.angle_table.add_perpendicular(rel)
    
    # checks sameclock
    def is_sameclock(self, p1: Point, p2: Point, p3: Point, p4: Point, p5: Point, p6: Point) -> bool:
        return ((p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x)) * \
                ((p5.x - p4.x) * (p6.y - p4.y) - (p5.y - p4.y) * (p6.x - p4.x)) > 0
    
    # checks if two segments are congruent via AR table
    def are_segments_congruent(self, seg1: frozenset, seg2: frozenset) -> Tuple[bool, List[RelationNode]]:
        if seg1 == seg2:
            return True, []
        if seg1 not in self.ratio_table.col_id or seg2 not in self.ratio_table.col_id:
            return False, []
        cong_row = [0] * self.ratio_table.table_length()
        cong_row[self.ratio_table.col_id[seg1]] = 1
        cong_row[self.ratio_table.col_id[seg2]] = -1
        
        is_spanned, parents = self.ratio_table.is_spanned(cong_row)
        return is_spanned, parents

    # checks if two angles are equal via AR table
    def are_angles_equal(self, angle1: Tuple[frozenset, frozenset], angle2: Tuple[frozenset, frozenset]) -> Tuple[bool, List[RelationNode]]:
        seg1_1, seg1_2 = angle1
        seg2_1, seg2_2 = angle2
        if not self.angles_look_equal(angle1, angle2):
            return False, []

        if (seg1_1 not in self.angle_table.col_id or
            seg1_2 not in self.angle_table.col_id or
            seg2_1 not in self.angle_table.col_id or
            seg2_2 not in self.angle_table.col_id):
            return False, []
        
        angle_row = [0] * self.angle_table.table_length()
        angle_row[self.angle_table.col_id[seg1_1]] += 1
        angle_row[self.angle_table.col_id[seg1_2]] += -1
        angle_row[self.angle_table.col_id[seg2_1]] += -1
        angle_row[self.angle_table.col_id[seg2_2]] += 1

        is_spanned, parents = self.angle_table.is_spanned(angle_row)
        return is_spanned, parents

    def are_equal_ratio(self, a: frozenset, b: frozenset, c: frozenset, d: frozenset) -> Tuple[bool, List[RelationNode]]:
        """Check if a/b = c/d via the RatioTable."""
        if (a not in self.ratio_table.col_id or b not in self.ratio_table.col_id or
            c not in self.ratio_table.col_id or d not in self.ratio_table.col_id):
            return False, []
        row = [0] * self.ratio_table.table_length()
        row[self.ratio_table.col_id[a]] += 1
        row[self.ratio_table.col_id[b]] -= 1
        row[self.ratio_table.col_id[c]] -= 1
        row[self.ratio_table.col_id[d]] += 1

        is_spanned, parents = self.ratio_table.is_spanned(row)
        return is_spanned, parents

    def are_points_collinear(self, p1: Point, p2: Point, p3: Point) -> Tuple[bool, List[RelationNode]]:
        """Check if three points are collinear via the AngleTable."""
        if (frozenset({p1, p2}) not in self.angle_table.col_id or
            frozenset({p2, p3}) not in self.angle_table.col_id or
            frozenset({p1, p3}) not in self.angle_table.col_id):
            return False, []

        gradient1 = (p2.y - p1.y) / (p2.x - p1.x)
        gradient2 = (p3.y - p2.y) / (p3.x - p2.x)
        if abs(gradient1 - gradient2) > 1e-6:
            return False, []

        row1 = [0] * self.angle_table.table_length()
        row1[self.angle_table.col_id[frozenset({p1, p2})]] += 1
        row1[self.angle_table.col_id[frozenset({p2, p3})]] += -1
        row2 = [0] * self.angle_table.table_length()
        row2[self.angle_table.col_id[frozenset({p1, p2})]] += 1
        row2[self.angle_table.col_id[frozenset({p1, p3})]] += -1
        
        is_spanned1, parents1 = self.angle_table.is_spanned(row1)
        is_spanned2, parents2 = self.angle_table.is_spanned(row2)
        if is_spanned1 and is_spanned2:
            return True, parents1 + parents2
        return False, []

    def are_points_parallel(self, p1: Point, p2: Point, p3: Point, p4: Point) -> Tuple[bool, List[RelationNode]]:
        """Check if two lines (p1,p2) and (p3,p4) are parallel via the AngleTable."""
        if (frozenset({p1, p2}) not in self.angle_table.col_id or
            frozenset({p3, p4}) not in self.angle_table.col_id):
            return False, []

        gradient1 = (p2.y - p1.y) / (p2.x - p1.x)
        gradient2 = (p4.y - p3.y) / (p4.x - p3.x)
        if abs(gradient1 - gradient2) > 1e-6:
            return False, []
        
        seg1 = frozenset({p1, p2})
        seg2 = frozenset({p3, p4})
        row = [0] * self.angle_table.table_length()
        row[self.angle_table.col_id[seg1]] = 1
        row[self.angle_table.col_id[seg2]] = -1

        is_spanned, parents = self.angle_table.is_spanned(row)
        return is_spanned, parents

    def are_points_perpendicular(self, p1: Point, p2: Point, p3: Point, p4: Point) -> Tuple[bool, List[RelationNode]]:
        """Check if two lines (p1,p2) and (p3,p4) are perpendicular via the AngleTable."""
        if (frozenset({p1, p2}) not in self.angle_table.col_id or
            frozenset({p3, p4}) not in self.angle_table.col_id):
            return False, []

        gradient1 = (p2.y - p1.y) / (p2.x - p1.x)
        gradient2 = (p4.y - p3.y) / (p4.x - p3.x)
        if abs(gradient1 * gradient2 + 1) > 1e-6:
            return False, []

        seg1 = frozenset({p1, p2})
        seg2 = frozenset({p3, p4})
        row1 = [0] * self.angle_table.table_length()
        row1[self.angle_table.col_id[seg1]] = 1
        row1[self.angle_table.col_id[seg2]] = -1  # Perpendicularity represented differently
        row1[self.angle_table.col_id[CONST_90]] = 1

        row2 = [0] * self.angle_table.table_length()
        row2[self.angle_table.col_id[seg1]] = 1
        row2[self.angle_table.col_id[seg2]] = -1  # Perpendicularity represented differently
        row2[self.angle_table.col_id[CONST_90]] = -1

        is_spanned1, parents1 = self.angle_table.is_spanned(row1)
        is_spanned2, parents2 = self.angle_table.is_spanned(row2)
        if is_spanned1 and is_spanned2:
            return True, parents1 + parents2
        elif is_spanned1:
            return True, parents1
        elif is_spanned2:
            return True, parents2
        return False, []

    def is_midpoint(self, mid: Point, p1: Point, p2: Point) -> Tuple[bool, List[RelationNode]]:
        """Check if mid is the midpoint of segment p1p2 via the RatioTable."""
        seg1 = frozenset({mid, p1})
        seg2 = frozenset({mid, p2})

        is_congruent, parents_cong = self.are_segments_congruent(seg1, seg2)
        is_collinear, parents_col = self.are_points_collinear(mid, p1, p2)
        if is_congruent and is_collinear:
            return True, parents_cong + parents_col
        return False, []

    def form_triangle(self, points: List[Point]) -> Optional[Tuple[Point, Point, Point, Point]]:
        """Given four points from a congruence relation, return them as triangle points; 
        if impossible, return None."""
        p1, p2, p3, p4 = points
        if p1 == p3:
            pass
        elif p1 == p4:
            p3, p4 = p4, p3
        elif p2 == p3:
            p1, p2 = p2, p1
        elif p2 == p4:
            p1, p2, p3, p4 = p2, p1, p4, p3
        else:
            return None, None, None, None
        
        return p1, p2, p3, p4  # p1 = p3
    
    def angles_look_equal(self, angle1: Tuple[frozenset, frozenset], angle2: Tuple[frozenset, frozenset]) -> bool:
        """Heuristic check if two angles look equal based on points."""
        seg1_1, seg1_2 = angle1
        seg2_1, seg2_2 = angle2
        p1, p2, p4, p3 = self.form_triangle(list(seg1_1) + list(seg1_2))
        q1, q2, q4, q3 = self.form_triangle(list(seg2_1) + list(seg2_2))

        # p1 is vertex of angle1, q1 is vertex of angle2
        # angle1 = p2 p1 p3, angle2 = q2 q1 q3
        angle1 = math.atan2(p2.y - p1.y, p2.x - p1.x) - math.atan2(p3.y - p1.y, p3.x - p1.x)
        angle2 = math.atan2(q2.y - q1.y, q2.x - q1.x) - math.atan2(q3.y - q1.y, q3.x - q1.x)
        angle1 = (angle1 + math.pi) % math.pi
        angle2 = (angle2 + math.pi) % math.pi
        return abs(angle1 - angle2) < 1e-6