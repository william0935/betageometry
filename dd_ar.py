import itertools
from relations import *
from problem import Problem
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
            # self.congABCD_congCDEF__congABEF,
            # self.eqangleABCDEF_eqangleDEFGHI__eqangleABCGHI,
            # self.eqangleABCDEF_eqangleFEDIHG__eqangleABCGHI,
            self.colABC__eqangleACDBCD_eqangleABDCBD_eqangleBADCAD,
            self.congMAMB_colMAB__midpointMAB,
            self.check_triangle_congruence_SSS,
            self.check_triangle_congruence_SAS,
            self.check_triangle_congruence_ASA,
            self.check_triangle_similarity_AA,
            self.check_triangle_similarity_SAS,
            self.check_triangle_similarity_SSS,
            self.eqratioABCDEFGH__eqratioABEFCDGH,
            self.congABEF_eqratioABCDEFGH__congCDGH,
            self.eqratioABACDEDF_colABC_colDEF__eqratioABBCDEEF
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
        for i, row in enumerate(table.rows):
            row_str = "  Row " + f"{i}: " + "".join(f"{val:8.2f}" if abs(val) > 1e-10 else "    0.00" for val in row)
            print(row_str)
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
        for col in self.problem.relations.get("col", []):
            self.angle_table.add_collinear(col)

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
                
        return self.problem.is_solved()
    
    def colABC__eqangleACDBCD_eqangleABDCBD_eqangleBADCAD(self) -> List[RelationNode]:
        new_relations = []
        collinears = self.problem.relations.get("col", [])
        
        for col_rel in collinears:
            if len(col_rel.points) != 3:
                continue
            for p4 in self.problem.points:
                if p4 not in col_rel.points:
                    p1, p2, p3 = col_rel.points
                    new_rel = [
                        EqualAngle(p1, p3, p4, p2, p3, p4, parents=[col_rel], rule="col_ABC__eqangle_ACDBCD"),
                        EqualAngle(p1, p2, p4, p3, p2, p4, parents=[col_rel], rule="col_ABC__eqangle_ABDCBD"),
                        EqualAngle(p2, p1, p4, p3, p1, p4, parents=[col_rel], rule="col_ABC__eqangle_BADCAD")
                    ]
                    new_relations.extend(new_rel)
                    
                    for eq_angle in new_rel:
                        self.angle_table.add_eqangle(eq_angle)

        return new_relations

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
            if seg_MA in self.angle_table.col_id and seg_MB in self.angle_table.col_id:
                col_row = [0] * self.angle_table.table_length()
                col_row[self.angle_table.col_id[seg_MA]] = 1
                col_row[self.angle_table.col_id[seg_MB]] = -1
                
                if self.angle_table.is_spanned(col_row):
                    col = next((c for c in collinears if set(c.points) == set({M, A, B})), None)
                    if col:
                        new_relations.append(Midpoint(
                            M, A, B,
                            parents=[cong, col],
                            rule="cong_MAMB_col_MAB__midpoint_MAB"
                        ))
        
        return new_relations
    
    def check_triangle_congruence_SSS(self) -> List[RelationNode]:
        """Check for triangle congruence using SSS criterion with caching."""
        new_relations = []
        triangles = list(itertools.combinations(self.problem.points, 3))
        
        for i in range(len(triangles)):
            for j in range(i + 1, len(triangles)):
                tri1 = triangles[i]
                tri2 = triangles[j]
                
                if set(tri1) == set(tri2):
                    continue
                
                cache_key = frozenset([frozenset(tri1), frozenset(tri2)])
                if cache_key in self.triangle_pair_cache and self.triangle_pair_cache[cache_key].get('SSS_checked'):
                    continue
                
                if cache_key not in self.triangle_pair_cache:
                    self.triangle_pair_cache[cache_key] = {}
                self.triangle_pair_cache[cache_key]['SSS_checked'] = True
                
                for perm in itertools.permutations(tri2):
                    p1, p2, p3 = tri1
                    p4, p5, p6 = perm
                    side1_1, side1_2, side1_3 = frozenset({p1, p2}), frozenset({p2, p3}), frozenset({p1, p3})
                    side2_1, side2_2, side2_3 = frozenset({p4, p5}), frozenset({p5, p6}), frozenset({p4, p6})
                    
                    if (self.are_segments_congruent(side1_1, side2_1) and
                        self.are_segments_congruent(side1_2, side2_2) and
                        self.are_segments_congruent(side1_3, side2_3)):
                        
                        sameclock = self.is_sameclock(p1, p2, p3, p4, p5, p6)
                        if sameclock:
                            new_relations.append(CongruentTriangle1(
                                p1, p2, p3, p4, p5, p6,
                                parents=[],
                                rule="check_triangle_congruence_SSS_sameclock"
                            ))
                        else:
                            new_relations.append(CongruentTriangle2(
                                p1, p2, p3, p4, p5, p6,
                                parents=[],
                                rule="check_triangle_congruence_SSS_opposite"
                            ))
                        break
        
        return new_relations
    
    def check_triangle_congruence_SAS(self) -> List[RelationNode]:
        """Check for triangle congruence using SAS criterion with caching."""
        new_relations = []
        triangles = list(itertools.combinations(self.problem.points, 3))
        
        for i in range(len(triangles)):
            for j in range(i + 1, len(triangles)):
                tri1 = triangles[i]
                tri2 = triangles[j]
                
                if set(tri1) == set(tri2):
                    continue
                
                cache_key = frozenset([frozenset(tri1), frozenset(tri2)])
                if cache_key in self.triangle_pair_cache and self.triangle_pair_cache[cache_key].get('SAS_checked'):
                    continue
                
                if cache_key not in self.triangle_pair_cache:
                    self.triangle_pair_cache[cache_key] = {}
                self.triangle_pair_cache[cache_key]['SAS_checked'] = True
                
                for perm in itertools.permutations(tri2):
                    p1, p2, p3 = tri1
                    p4, p5, p6 = perm
                    side1_1, side1_2 = frozenset({p1, p2}), frozenset({p2, p3})
                    side2_1, side2_2 = frozenset({p4, p5}), frozenset({p5, p6})
                    angle1 = (frozenset({p1, p2}), frozenset({p2, p3}))
                    angle2 = (frozenset({p4, p5}), frozenset({p5, p6}))
                    
                    if (self.are_segments_congruent(side1_1, side2_1) and
                        self.are_segments_congruent(side1_2, side2_2) and
                        self.are_angles_equal(angle1, angle2)):
                        
                        sameclock = self.is_sameclock(p1, p2, p3, p4, p5, p6)
                        if sameclock:
                            new_relations.append(CongruentTriangle1(
                                p1, p2, p3, p4, p5, p6,
                                parents=[],
                                rule="check_triangle_congruence_SAS_sameclock"
                            ))
                        else:
                            new_relations.append(CongruentTriangle2(
                                p1, p2, p3, p4, p5, p6,
                                parents=[],
                                rule="check_triangle_congruence_SAS_opposite"
                            ))
                        break
        
        return new_relations
    
    def check_triangle_congruence_ASA(self) -> List[RelationNode]:
        """Check for triangle congruence using ASA/AAS criterion with caching."""
        new_relations = []
        triangles = list(itertools.combinations(self.problem.points, 3))
        
        for i in range(len(triangles)):
            for j in range(i + 1, len(triangles)):
                tri1 = triangles[i]
                tri2 = triangles[j]
                
                if set(tri1) == set(tri2):
                    continue
                
                cache_key = frozenset([frozenset(tri1), frozenset(tri2)])
                if cache_key in self.triangle_pair_cache and self.triangle_pair_cache[cache_key].get('ASA_checked'):
                    continue
                
                if cache_key not in self.triangle_pair_cache:
                    self.triangle_pair_cache[cache_key] = {}
                self.triangle_pair_cache[cache_key]['ASA_checked'] = True
                
                for perm in itertools.permutations(tri2):
                    p1, p2, p3 = tri1
                    p4, p5, p6 = perm
                    angles1 = [
                        (frozenset({p1, p2}), frozenset({p2, p3})),
                        (frozenset({p2, p3}), frozenset({p3, p1})),
                        (frozenset({p3, p1}), frozenset({p1, p2}))
                    ]
                    angles2 = [
                        (frozenset({p4, p5}), frozenset({p5, p6})),
                        (frozenset({p5, p6}), frozenset({p6, p4})),
                        (frozenset({p6, p4}), frozenset({p4, p5}))
                    ]
                    
                    equal_angles = 0
                    for a1, a2 in zip(angles1, angles2):
                        if self.are_angles_equal(a1, a2):
                            equal_angles += 1
                    
                    if equal_angles >= 2:
                        sides1 = [frozenset({p1, p2}), frozenset({p2, p3}), frozenset({p1, p3})]
                        sides2 = [frozenset({p4, p5}), frozenset({p5, p6}), frozenset({p4, p6})]
                        has_one_side = any(self.are_segments_congruent(s1, s2) 
                                         for s1 in sides1 for s2 in sides2)
                        
                        if has_one_side:
                            sameclock = self.is_sameclock(p1, p2, p3, p4, p5, p6)
                            if sameclock:
                                new_relations.append(CongruentTriangle1(
                                    p1, p2, p3, p4, p5, p6,
                                    parents=[],
                                    rule="check_triangle_congruence_ASA_sameclock"
                                ))
                            else:
                                new_relations.append(CongruentTriangle2(
                                    p1, p2, p3, p4, p5, p6,
                                    parents=[],
                                    rule="check_triangle_congruence_ASA_opposite"
                                ))
                            break
        
        return new_relations
    
    def check_triangle_similarity_AA(self) -> List[RelationNode]:
        """Check for triangle similarity using AA criterion with caching."""
        new_relations = []
        triangles = list(itertools.combinations(self.problem.points, 3))
        
        for i in range(len(triangles)):
            for j in range(i + 1, len(triangles)):
                tri1 = triangles[i]
                tri2 = triangles[j]
                
                if set(tri1) == set(tri2):
                    continue
                
                cache_key = frozenset([frozenset(tri1), frozenset(tri2)])
                if cache_key in self.triangle_pair_cache and self.triangle_pair_cache[cache_key].get('AA_checked'):
                    continue
                
                if cache_key not in self.triangle_pair_cache:
                    self.triangle_pair_cache[cache_key] = {}
                self.triangle_pair_cache[cache_key]['AA_checked'] = True
                
                for perm in itertools.permutations(tri2):
                    p1, p2, p3 = tri1
                    p4, p5, p6 = perm
                    
                    angles1 = [
                        (frozenset({p1, p2}), frozenset({p2, p3})),
                        (frozenset({p2, p3}), frozenset({p3, p1})),
                        (frozenset({p3, p1}), frozenset({p1, p2}))
                    ]
                    angles2 = [
                        (frozenset({p4, p5}), frozenset({p5, p6})),
                        (frozenset({p5, p6}), frozenset({p6, p4})),
                        (frozenset({p6, p4}), frozenset({p4, p5}))
                    ]
                    
                    equal_angles = 0
                    for a1, a2 in zip(angles1, angles2):
                        if self.are_angles_equal(a1, a2):
                            equal_angles += 1
                    
                    if equal_angles >= 2:
                        sameclock = self.is_sameclock(p1, p2, p3, p4, p5, p6)
                        if sameclock:
                            new_relations.append(SimilarTriangle1(
                                p1, p2, p3, p4, p5, p6,
                                parents=[],
                                rule="check_triangle_similarity_AA_sameclock"
                            ))
                        else:
                            new_relations.append(SimilarTriangle2(
                                p1, p2, p3, p4, p5, p6,
                                parents=[],
                                rule="check_triangle_similarity_AA_opposite"
                            ))
                        break
        
        return new_relations
    
    def check_triangle_similarity_SAS(self) -> List[RelationNode]:
        """Check for triangle similarity using SAS criterion with caching."""
        new_relations = []
        triangles = list(itertools.combinations(self.problem.points, 3))
        
        for i in range(len(triangles)):
            for j in range(i + 1, len(triangles)):
                tri1 = triangles[i]
                tri2 = triangles[j]
                
                if set(tri1) == set(tri2):
                    continue
                
                cache_key = frozenset([frozenset(tri1), frozenset(tri2)])
                if cache_key in self.triangle_pair_cache and self.triangle_pair_cache[cache_key].get('SAS_sim_checked'):
                    continue
                
                if cache_key not in self.triangle_pair_cache:
                    self.triangle_pair_cache[cache_key] = {}
                self.triangle_pair_cache[cache_key]['SAS_sim_checked'] = True
                
                for perm in itertools.permutations(tri2):
                    p1, p2, p3 = tri1
                    p4, p5, p6 = perm
                    angle1 = (frozenset({p1, p2}), frozenset({p2, p3}))
                    angle2 = (frozenset({p4, p5}), frozenset({p5, p6}))
                    
                    if self.are_angles_equal(angle1, angle2):
                        # Check if sides are proportional: p1p2/p4p5 = p2p3/p5p6
                        # This means p1p2 * p5p6 = p2p3 * p4p5 in ratio table
                        # For now, skip this complex check - would need eqratio in AR table
                        pass
        
        return new_relations
    
    def check_triangle_similarity_SSS(self) -> List[RelationNode]:
        """Check for triangle similarity using SSS criterion with caching."""
        new_relations = []
        # This requires checking if all three side ratios are equal
        # Would need more sophisticated ratio checks in AR table
        # Skip for now as it's complex to implement
        return new_relations
    
    def eqratioABCDEFGH__eqratioABEFCDGH(self) -> List[RelationNode]:
        new_relations = []
        eqratios = self.problem.relations.get("eqratio", [])
        
        for eqr in eqratios:
            p1, p2, p3, p4, p5, p6, p7, p8 = eqr.points
            new_relations.append(EqualRatio(
                p1, p2, p5, p6, p3, p4, p7, p8,
                parents=[eqr],
                rule="eqratioABCDEFGH__eqratioABEFCDGH"
            ))

        return new_relations
    
    def congABEF_eqratioABCDEFGH__congCDGH(self) -> List[RelationNode]:
        new_relations = []
        eqratios = self.problem.relations.get("eqratio", [])
        
        for eqr in eqratios:
            p1, p2, p3, p4, p5, p6, p7, p8 = eqr.points
            segment_pairs = [
                (frozenset({p1, p2}), frozenset({p5, p6}), frozenset({p3, p4}), frozenset({p7, p8})),
                (frozenset({p3, p4}), frozenset({p7, p8}), frozenset({p1, p2}), frozenset({p5, p6})),
                (frozenset({p1, p2}), frozenset({p3, p4}), frozenset({p5, p6}), frozenset({p7, p8})),
                (frozenset({p5, p6}), frozenset({p7, p8}), frozenset({p1, p2}), frozenset({p3, p4}))
            ]
            
            for seg1, seg2, seg3, seg4 in segment_pairs:
                if self.are_segments_congruent(seg1, seg2):
                    new_relations.append(Congruent(
                        list(seg3)[0], list(seg3)[1], list(seg4)[0], list(seg4)[1],
                        parents=[eqr],
                        rule="congABEF_eqratioABCDEFGH__congCDGH"
                    ))
                    break

        return new_relations
    
    def paraABCD_colACE_colBDE__eqratioACCEBDDE(self) -> List[RelationNode]:
        new_relations = []
        paralles = self.problem.relations.get("para", [])
        collinears = self.problem.relations.get("col", [])
        for para in paralles:
            p1, p2, p3, p4 = para.points
            for p5 in self.problem.points:
                if p5 in {p1, p2, p3, p4}:
                    continue
                col1 = None
                col2 = None
                col3 = None
                col4 = None
                for col in collinears:
                    if set({p1, p3, p5}) == set(col.points):
                        col1 = col
                    if set({p2, p4, p5}) == set(col.points):
                        col2 = col
                    if set({p1, p4, p5}) == set(col.points):
                        col3 = col
                    if set({p2, p3, p5}) == set(col.points):
                        col4 = col
                if col3 and col4:
                    p3, p4 = p4, p3
                elif col1 and col2:
                    pass
                else:
                    continue

                new_relations.append(EqualRatio(
                    p1, p3, p3, p5, p2, p4, p4, p5,
                    parents=[para, col1, col2],
                    rule="paraABCD_colACE_colBDE__eqratioACCEBDDE"
                ))

        return new_relations
    
    def eqratioABACDEDF_colABC_colDEF__eqratioABBCDEEF(self) -> List[RelationNode]:
        new_relations = []
        eqratios = self.problem.relations.get("eqratio", [])
        collinears = self.problem.relations.get("col", [])
        for eqr in eqratios:
            p1, p2, p3, p4, p5, p6, p7, p8 = eqr.points
            if p1 == p3:
                pass
            elif p1 == p4:
                p3, p4 = p4, p3
            elif p2 == p3:
                p1, p2 = p2, p1
            elif p2 == p4:
                p1, p2, p3, p4 = p2, p1, p4, p3
            else:
                continue

            if p5 == p7:
                pass
            elif p5 == p8:
                p7, p8 = p8, p7
            elif p6 == p7:
                p5, p6 = p6, p5
            elif p6 == p8:
                p5, p6, p7, p8 = p6, p5, p8, p7
            else:
                continue

            
            col1 = None
            col2 = None
            for col in collinears:
                if set({p1, p2, p4}) == set(col.points):
                    col1 = col
                if set({p5, p6, p8}) == set(col.points):
                    col2 = col
            
            if not col1 or not col2:
                continue

            new_relations.extend([EqualRatio(
                p1, p2, p2, p4, p5, p6, p6, p8,
                parents=[eqr, col1, col2],
                rule="eqratioABACDEDF_colABC_colDEF__eqratioABBCDEEF"
            ), EqualRatio(
                p1, p4, p2, p4, p5, p8, p6, p8,
                parents=[eqr, col1, col2],
                rule="eqratioABACDEDF_colABC_colDEF__eqratioABBCDEEF"
            )])

        return new_relations



    """
    Helper methods for specific rules starts here.
    """
    
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
    
    # checks sameclock
    def is_sameclock(self, p1: Point, p2: Point, p3: Point, p4: Point, p5: Point, p6: Point) -> bool:
        return ((p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x)) * \
                ((p5.x - p4.x) * (p6.y - p4.y) - (p5.y - p4.y) * (p6.x - p4.x)) > 0
    
    # checks if two segments are congruent via AR table
    def are_segments_congruent(self, seg1: frozenset, seg2: frozenset) -> bool:
        if seg1 == seg2:
            return True
        if seg1 not in self.ratio_table.col_id or seg2 not in self.ratio_table.col_id:
            return False
        cong_row = [0] * self.ratio_table.table_length()
        cong_row[self.ratio_table.col_id[seg1]] = 1
        cong_row[self.ratio_table.col_id[seg2]] = -1
        return self.ratio_table.is_spanned(cong_row)
    
    # checks if two angles are equal via AR table
    def are_angles_equal(self, angle1: Tuple[frozenset, frozenset], angle2: Tuple[frozenset, frozenset]) -> bool:
        seg1_1, seg1_2 = angle1
        seg2_1, seg2_2 = angle2
        if (seg1_1 not in self.angle_table.col_id or seg1_2 not in self.angle_table.col_id or
            seg2_1 not in self.angle_table.col_id or seg2_2 not in self.angle_table.col_id):
            return False
        angle_row = [0] * self.angle_table.table_length()
        angle_row[self.angle_table.col_id[seg1_1]] += 1
        angle_row[self.angle_table.col_id[seg1_2]] += -1
        angle_row[self.angle_table.col_id[seg2_1]] += -1
        angle_row[self.angle_table.col_id[seg2_2]] += 1
        return self.angle_table.is_spanned(angle_row)