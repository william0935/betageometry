import itertools
from relations import *
from Problem import Problem
from typing import List, Tuple, Optional


class DeductiveDatabase:    
    def __init__(self, problem: Problem):
        self.problem = problem

        # Rule database - each rule is a method that returns List[RelationNode]
        # ADD THE NEW RULES HERE
        self.rules = [
            self.congABCD_congCDEF__congABEF,
            self.eqangleABCDEF_eqangleDEFGHI__eqangleABCGHI,
            self.eqangleABCDEF_eqangleFEDIHG__eqangleABCGHI,
            self.colABC__eqangleACDBCD_eqangleABDCBD_eqangleBADCAD,
            self.congMAMB_colMAB__midpMAB,
            self.congABDE_congBCEF_eqangleABCDEF_sameclockABCDEF__contri1ABCDEF,
            self.congABDE_congBCEF_eqangleABCFED_sameclockABCFED__contri2ABCDEF,
            self.congABDE_eqangleABCDEF_eqangleCABFDE_sameclockABCDEF__contri1ABCDEF,
            self.congABDE_eqangleABCFED_eqangleCABDFE_sameclockABCFED__contri2ABCDEF,
            self.congABDE_congBCEF_congCAFD__contriABCDEF,
            self.eqangleABCDEF_eqangleCABFDE_sameclockABCDEF__simtri1ABCDEF,
            self.eqangleABCFED_eqangleCABDFE_sameclockABCFED__simtri2ABCDEF,
            self.eqratioABCDEFGH__eqratioABEFCDGH,
            self.congABEF_eqratioABCDEFGH__congCDGH,
            self.paraABCD_colACE_colBDE__eqratioACCEBDDE,
            self.eqratioABACDEDF_colABC_colDEF__eqratioABBCDEEF,
            self.congOAOB_congOBOC_colOAB__perpACBC,
            self.congOAOB_congOAOC_congOAOD__cyclicABCD,
            self.midpMAB_midpNAC__paraMNBC
        ]

    
    def apply_deduction_rules(self, max_iterations: int) -> bool:   
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
            
            # Check if we made progress
            if not progress_made:
                print(f"No new relations in iteration {iteration}. Stopping.")
                break  # No new relations derived, stop
                
            if self.problem.is_solved():
                return True
                
        return self.problem.is_solved()
    
    # Rule: transitivity of congruent segments
    def congABCD_congCDEF__congABEF(self) -> List[RelationNode]:
        new_relations = []
        congruences = self.problem.relations.get("cong", [])
        
        for i in range(len(congruences)):
            for j in range(i + 1, len(congruences)):
                p1, p2, p3, p4 = congruences[i].points
                p5, p6, p7, p8 = congruences[j].points
                seg1 = frozenset({p1, p2})
                seg2 = frozenset({p3, p4})
                seg3 = frozenset({p5, p6})
                seg4 = frozenset({p7, p8})
                if seg1 == seg4:
                    p5, p6, p7, p8 = p7, p8, p5, p6
                elif seg2 == seg3:
                    p1, p2, p3, p4 = p3, p4, p1, p2
                elif seg2 == seg4:
                    p1, p2, p3, p4, p5, p6, p7, p8 = p3, p4, p1, p2, p7, p8, p5, p6
                elif seg1 == seg3:
                    pass
                else:
                    continue
                
                if frozenset({p3, p4}) == frozenset({p7, p8}):
                    continue

                new_relations.append(Congruent(
                    p3, p4, p7, p8, 
                    parents=[congruences[i], congruences[j]], 
                    rule="congABCD_congCDEF__congABEF"
                ))

        return new_relations
    
    # Rule: transitivity of equal angles
    def eqangleABCDEF_eqangleDEFGHI__eqangleABCGHI(self) -> List[RelationNode]:
        new_relations = []
        eqangles = self.problem.relations.get("eqangle", [])
        
        for i in range(len(eqangles)):
            for j in range(i + 1, len(eqangles)):
                p1, p2, p3, p4, p5, p6 = eqangles[i].points
                p7, p8, p9, p10, p11, p12 = eqangles[j].points
                angle1 = (frozenset({p1, p2}), frozenset({p2, p3}))
                angle2 = (frozenset({p4, p5}), frozenset({p5, p6}))
                angle3 = (frozenset({p7, p8}), frozenset({p8, p9}))
                angle4 = (frozenset({p10, p11}), frozenset({p11, p12}))
                if angle1 == angle4:
                    p7, p8, p9, p10, p11, p12 = p10, p11, p12, p7, p8, p9
                elif angle2 == angle3:
                    p1, p2, p3, p4, p5, p6 = p4, p5, p6, p1, p2, p3
                elif angle2 == angle4:
                    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12 = p4, p5, p6, p1, p2, p3, p10, p11, p12, p7, p8, p9
                elif angle1 == angle3:
                    pass
                else:
                    continue
                
                if ({frozenset({p4, p5}), frozenset({p5, p6})}) == ({frozenset({p10, p11}), frozenset({p11, p12})}):
                    continue

                new_relations.append(EqualAngle(
                    p4, p5, p6, p10, p11, p12,
                    parents=[eqangles[i], eqangles[j]], 
                    rule="eqangleABCDEF_eqangleDEFGHI__eqangleABCGHI"
                ))

        return new_relations
    
    # Rule: transitivity of equal angles (reversed second angle)
    def eqangleABCDEF_eqangleFEDIHG__eqangleABCGHI(self) -> List[RelationNode]:
        new_relations = []
        eqangles = self.problem.relations.get("eqangle", [])

        for i in range(len(eqangles)):
            for j in range(i + 1, len(eqangles)):
                p1, p2, p3, p4, p5, p6 = eqangles[i].points
                p9, p8, p7, p12, p11, p10 = eqangles[j].points
                angle1 = (frozenset({p1, p2}), frozenset({p2, p3}))
                angle2 = (frozenset({p4, p5}), frozenset({p5, p6}))
                angle3 = (frozenset({p7, p8}), frozenset({p8, p9}))
                angle4 = (frozenset({p10, p11}), frozenset({p11, p12}))
                
                if angle1 == angle4:
                    p7, p8, p9, p10, p11, p12 = p10, p11, p12, p7, p8, p9
                elif angle2 == angle3:
                    p1, p2, p3, p4, p5, p6 = p4, p5, p6, p1, p2, p3
                elif angle2 == angle4:
                    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12 = p4, p5, p6, p1, p2, p3, p10, p11, p12, p7, p8, p9
                elif angle1 == angle3:
                    pass
                else:
                    continue

                if ({frozenset({p4, p5}), frozenset({p5, p6})}) == ({frozenset({p10, p11}), frozenset({p11, p12})}):
                    continue

                new_relations.append(EqualAngle(
                    p4, p5, p6, p10, p11, p12,
                    parents=[eqangles[i], eqangles[j]], 
                    rule="eqangleABCDEF_eqangleFEDIHG__eqangleABCGHI"
                ))
        
        return new_relations
    
    # Rule: collinear points -> equal angles
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
                        EqualAngle(p1, p3, p4, p2, p3, p4, parents=[col_rel], rule="colABC__eqangleACDBCD"),
                        EqualAngle(p1, p2, p4, p3, p2, p4, parents=[col_rel], rule="colABC__eqangleABDCBD"),
                        EqualAngle(p2, p1, p4, p3, p1, p4, parents=[col_rel], rule="colABC__eqangleBADCAD")
                    ]
                    new_relations.extend(new_rel)

        return new_relations

    # Rule: congruent segments + collinear -> midpoint
    def congMAMB_colMAB__midpMAB(self) -> List[RelationNode]:
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

            for col in collinears:
                if set({p1, p2, p4}) == set(col.points):
                    new_relations.append(Midpoint(
                        p1, p2, p4,
                        parents=[cong, col],
                        rule="congMAMB_colMAB__midpMAB"
                    ))
        
        return new_relations
    
    # Rule: triangle congruence SAS with same orientation
    def congABDE_congBCEF_eqangleABCDEF_sameclockABCDEF__contri1ABCDEF(self) -> List[RelationNode]:
        new_relations = []
        eqangles = self.problem.relations.get("eqangle", [])
        
        for eq in eqangles:
            p1, p2, p3, p4, p5, p6 = eq.points
            sameclock = True if self.is_sameclock(p1, p2, p3, p4, p5, p6) else False
            if not sameclock or (p1, p2, p3) == (p4, p5, p6):
                continue

            cong1_exists, cong1 = self.is_cong(p1, p2, p4, p5)
            cong2_exists, cong2 = self.is_cong(p2, p3, p5, p6)
            if cong1_exists and cong2_exists:
                parents = [eq]
                if cong1:
                    parents.append(cong1)
                if cong2:
                    parents.append(cong2)

                new_relations.append(CongruentTriangle1(
                    p1, p2, p3, p4, p5, p6, 
                    parents=parents,
                    rule="congABDE_congBCEF_eqangleABCDEF_sameclockABCDEF__contri1ABCDEF"
                ))

        return new_relations

    # Rule: triangle congruence SAS with opposite orientation
    def congABDE_congBCEF_eqangleABCFED_sameclockABCFED__contri2ABCDEF(self) -> List[RelationNode]:
        new_relations = []
        eqangles = self.problem.relations.get("eqangle", [])
        
        for eq in eqangles:
            p1, p2, p3, p6, p5, p4 = eq.points
            sameclock = True if self.is_sameclock(p1, p2, p3, p4, p5, p6) else False
            if sameclock or (p1, p2, p3) == (p4, p5, p6):
                continue

            cong1_exists, cong1 = self.is_cong(p1, p2, p4, p5)
            cong2_exists, cong2 = self.is_cong(p2, p3, p5, p6)
            if cong1_exists and cong2_exists:
                parents = [eq]
                if cong1:
                    parents.append(cong1)
                if cong2:
                    parents.append(cong2)

                new_relations.append(CongruentTriangle2(
                    p1, p2, p3, p4, p5, p6, 
                    parents=parents,
                    rule="congABDE_congBCEF_eqangleABCFED_sameclockABCFED__contri2ABCDEF"
                ))

        return new_relations
    
    # Rule: triangle Congruence ASA/AAS
    def congABDE_eqangleABCDEF_eqangleCABFDE_sameclockABCDEF__contri1ABCDEF(self) -> List[RelationNode]:
        new_relations = []
        eqangles = self.problem.relations.get("eqangle", [])
        
        for i in range(len(eqangles)):
            for j in range(i + 1, len(eqangles)):
                p1, p2, p3, p4, p5, p6 = eqangles[i].points
                p7, p8, p9, p10, p11, p12 = eqangles[j].points
                
                if {p1, p2, p3} == {p7, p8, p9} and {p4, p5, p6} == {p10, p11, p12}:
                    pass
                elif {p1, p2, p3} == {p10, p11, p12} and {p4, p5, p6} == {p7, p8, p9}:
                    p7, p8, p9, p10, p11, p12 = p10, p11, p12, p7, p8, p9
                else:
                    continue

                if frozenset({(p1, p4), (p2, p5), (p3, p6)}) != \
                   frozenset({(p7, p10), (p8, p11), (p9, p12)}) or \
                   p2 == p8 or p5 == p11 or \
                   not self.is_sameclock(p1, p2, p3, p4, p5, p6):
                    continue

                cong1_exists, cong1 = self.is_cong(p1, p2, p4, p5)
                cong2_exists, cong2 = self.is_cong(p2, p3, p5, p6)
                cong3_exists, cong3 = self.is_cong(p1, p3, p4, p6)
                parents = [eqangles[i], eqangles[j]]
                if cong1_exists:
                    if cong1:
                        parents.append(cong1)
                elif cong2_exists:
                    if cong2:
                        parents.append(cong2)
                elif cong3_exists:
                    if cong3:
                        parents.append(cong3)
                else:
                    continue

                new_relations.append(CongruentTriangle1(
                    p1, p2, p3, p4, p5, p6,
                    parents=parents,
                    rule="congABDE_eqangleABCDEF_eqangleCABFDE__contri1ABCDEF"
                ))

        return new_relations
    
    # Rule: triangle congruence ASA/AAS with opposite orientation
    def congABDE_eqangleABCFED_eqangleCABDFE_sameclockABCFED__contri2ABCDEF(self) -> List[RelationNode]:
        new_relations = []
        eqangles = self.problem.relations.get("eqangle", [])
        
        for i in range(len(eqangles)):
            for j in range(i + 1, len(eqangles)):
                p1, p2, p3, p6, p5, p4 = eqangles[i].points
                p7, p8, p9, p12, p11, p10 = eqangles[j].points
                
                if {p1, p2, p3} == {p7, p8, p9} and {p4, p5, p6} == {p10, p11, p12}:
                    pass
                elif {p1, p2, p3} == {p10, p11, p12} and {p4, p5, p6} == {p7, p8, p9}:
                    p7, p8, p9, p10, p11, p12 = p10, p11, p12, p7, p8, p9
                else:
                    continue

                if frozenset({(p1, p4), (p2, p5), (p3, p6)}) != \
                   frozenset({(p7, p10), (p8, p11), (p9, p12)}) or \
                   p2 == p8 or p5 == p11 or \
                   self.is_sameclock(p1, p2, p3, p4, p5, p6):
                    continue

                cong1_exists, cong1 = self.is_cong(p1, p2, p4, p5)
                cong2_exists, cong2 = self.is_cong(p2, p3, p5, p6)
                cong3_exists, cong3 = self.is_cong(p1, p3, p4, p6)
                parents = [eqangles[i], eqangles[j]]
                if cong1_exists:
                    if cong1:
                        parents.append(cong1)
                elif cong2_exists:
                    if cong2:
                        parents.append(cong2)
                elif cong3_exists:
                    if cong3:
                        parents.append(cong3)
                else:
                    continue

                new_relations.append(CongruentTriangle2(
                    p1, p2, p3, p4, p5, p6,
                    parents=parents,
                    rule="congABDE_eqangleABCFED_eqangleCBAFED__contri2ABCDEF"
                ))

        return new_relations
    
    # Rule: triangle congruence SSS
    def congABDE_congBCEF_congCAFD__contriABCDEF(self):
        new_relations = []
        congruences = self.problem.relations.get("cong", [])
        for i in range(len(congruences)):
            p1, p2, p3, p4 = congruences[i].points
            for j in range(i + 1, len(congruences)):
                p5, p6, p7, p8 = congruences[j].points
                if len({p1, p2, p5, p6}) == 3 and len({p3, p4, p7, p8}) == 3:
                    if p1 == p5:
                        pass
                    elif p1 == p6:
                        p5, p6 = p6, p5
                    elif p2 == p5:
                        p1, p2 = p2, p1
                    elif p2 == p6:
                        p1, p2, p5, p6 = p2, p1, p6, p5
                    else:
                        continue

                    if p3 == p7:
                        pass
                    elif p3 == p8:
                        p7, p8 = p8, p7
                    elif p4 == p7:
                        p3, p4 = p4, p3
                    elif p4 == p8:
                        p3, p4, p7, p8 = p4, p3, p8, p7
                    else:
                        continue

                    cong3_exists, cong3 = self.is_cong(p2, p6, p4, p8)
                    if cong3_exists:
                        parents = [congruences[i], congruences[j]]
                        if cong3:
                            parents.append(cong3)

                        if self.is_sameclock(p1, p2, p6, p3, p4, p8):
                            new_relations.append(CongruentTriangle1(
                                p1, p2, p6, p3, p4, p8,
                                parents=parents,
                                rule="congABDE_congBCEF_congCAFD_sameclockABCDEF__contri1ABCDEF"
                            ))
                        else:
                            new_relations.append(CongruentTriangle2(
                                p1, p2, p6, p3, p4, p8,
                                parents=parents,
                                rule="congABDE_congBCEF_congCAFD_sameclockABCFED__contri2ABCDEF"
                            ))

                p5, p6, p7, p8 = p7, p8, p5, p6
                if len({p1, p2, p5, p6}) == 3 and len({p3, p4, p7, p8}) == 3:
                    if p1 == p5:
                        pass
                    elif p1 == p6:
                        p5, p6 = p6, p5
                    elif p2 == p5:
                        p1, p2 = p2, p1
                    elif p2 == p6:
                        p1, p2, p5, p6 = p2, p1, p6, p5
                    else:
                        continue

                    if p3 == p7:
                        pass
                    elif p3 == p8:
                        p7, p8 = p8, p7
                    elif p4 == p7:
                        p3, p4 = p4, p3
                    elif p4 == p8:
                        p3, p4, p7, p8 = p4, p3, p8, p7
                    else:
                        continue

                    cong3_exists, cong3 = self.is_cong(p2, p6, p4, p8)
                    if cong3_exists:
                        parents = [congruences[i], congruences[j]]
                        if cong3:
                            parents.append(cong3)
                        
                        if self.is_sameclock(p1, p2, p6, p3, p4, p8):
                            new_relations.append(CongruentTriangle1(
                                p1, p2, p6, p3, p4, p8,
                                parents=parents,
                                rule="congABDE_congBCEF_congCAFD_sameclockABCDEF__contri1ABCDEF"
                            ))
                        else:
                            new_relations.append(CongruentTriangle2(
                                p1, p2, p6, p3, p4, p8,
                                parents=parents,
                                rule="congABDE_congBCEF_congCAFD_sameclockABCFED__contri2ABCDEF"
                            ))

        return new_relations

    # Rule: triangle similarity AA
    def eqangleABCDEF_eqangleCABFDE_sameclockABCDEF__simtri1ABCDEF(self) -> List[RelationNode]:
        new_relations = []
        eqangles = self.problem.relations.get("eqangle", [])
        
        for i in range(len(eqangles)):
            for j in range(i + 1, len(eqangles)):
                p1, p2, p3, p4, p5, p6 = eqangles[i].points
                p7, p8, p9, p10, p11, p12 = eqangles[j].points
                
                if {p1, p2, p3} == {p7, p8, p9} and {p4, p5, p6} == {p10, p11, p12}:
                    pass
                elif {p1, p2, p3} == {p10, p11, p12} and {p4, p5, p6} == {p7, p8, p9}:
                    p7, p8, p9, p10, p11, p12 = p10, p11, p12, p7, p8, p9
                else:
                    continue

                if frozenset({(p1, p4), (p2, p5), (p3, p6)}) != \
                   frozenset({(p7, p10), (p8, p11), (p9, p12)}) or \
                   p2 == p8 or p5 == p11 or \
                   not self.is_sameclock(p1, p2, p3, p4, p5, p6):
                    continue

                new_relations.append(SimilarTriangle1(
                    p1, p2, p3, p4, p5, p6,
                    parents=[eqangles[i], eqangles[j]],
                    rule="eqangleABCDEF_eqangleCABFDE_sameclockABCDEF__simtri1ABCDEF"
                ))

        return new_relations
    
    # Rule: triangle similarity AA with opposite orientation
    def eqangleABCFED_eqangleCABDFE_sameclockABCFED__simtri2ABCDEF(self) -> List[RelationNode]:
        new_relations = []
        eqangles = self.problem.relations.get("eqangle", [])
        
        for i in range(len(eqangles)):
            for j in range(i + 1, len(eqangles)):
                p1, p2, p3, p6, p5, p4 = eqangles[i].points
                p7, p8, p9, p12, p11, p10 = eqangles[j].points
                
                if {p1, p2, p3} == {p7, p8, p9} and {p4, p5, p6} == {p10, p11, p12}:
                    pass
                elif {p1, p2, p3} == {p10, p11, p12} and {p4, p5, p6} == {p7, p8, p9}:
                    p7, p8, p9, p10, p11, p12 = p10, p11, p12, p7, p8, p9
                else:
                    continue

                if frozenset({(p1, p4), (p2, p5), (p3, p6)}) != \
                   frozenset({(p7, p10), (p8, p11), (p9, p12)}) or \
                   p2 == p8 or p5 == p11 or \
                   self.is_sameclock(p1, p2, p3, p4, p5, p6):
                    continue

                new_relations.append(SimilarTriangle2(
                    p1, p2, p3, p4, p5, p6,
                    parents=[eqangles[i], eqangles[j]],
                    rule="eqangleABCFED_eqangleCABDFE_sameclockABCFED__simtri2ABCDEF"
                ))

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
    
    # Rule: Congruent segments + equal ratios -> congruent segments
    def congABEF_eqratioABCDEFGH__congCDGH(self) -> List[RelationNode]:
        new_relations = []
        congruences = self.problem.relations.get("cong", [])
        eqratios = self.problem.relations.get("eqratio", [])
        
        for eqr in eqratios:
            p1, p2, p3, p4, p5, p6, p7, p8 = eqr.points
            for cong in congruences:
                p9, p10, p11, p12 = cong.points
                if frozenset({frozenset({p1, p2}), frozenset({p5, p6})}) == \
                   frozenset({frozenset({p9, p10}), frozenset({p11, p12})}):
                    pass
                elif frozenset({frozenset({p3, p4}), frozenset({p7, p8})}) == \
                     frozenset({frozenset({p9, p10}), frozenset({p11, p12})}):
                    p1, p2, p3, p4, p5, p6, p7, p8 = p3, p4, p7, p8, p1, p2, p5, p6
                elif frozenset({frozenset({p1, p2}), frozenset({p3, p4})}) == \
                     frozenset({frozenset({p9, p10}), frozenset({p11, p12})}):
                    p3, p4, p5, p6 = p5, p6, p3, p4
                elif frozenset({frozenset({p5, p6}), frozenset({p7, p8})}) == \
                     frozenset({frozenset({p9, p10}), frozenset({p11, p12})}):
                    p1, p2, p5, p6 = p5, p6, p1, p2
                else:
                    continue

                new_relations.append(Congruent(
                    p3, p4, p7, p8,
                    parents=[eqr, cong],
                    rule="congABEF_eqratioABCDEFGH__congCDGH"
                ))

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
        
    def congOAOB_congOBOC_colOAB__perpACBC(self) -> List[RelationNode]:
        new_relations = []
        congruences = self.problem.relations.get("cong", [])
        collinears = self.problem.relations.get("col", [])
        
        for i in range(len(congruences)):
            cong1 = congruences[i]
            p1, p2, p3, p4 = cong1.points
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

            for j in range(i + 1, len(congruences)):
                cong2 = congruences[j]
                p5, p6, p7, p8 = cong2.points
                if len({p1, p2, p3, p4, p5, p6, p7, p8}) != 4:
                    continue

                if p5 == p8:
                    p7, p8 = p8, p7
                elif p6 == p7:
                    p5, p6 = p6, p5
                elif p6 == p8:
                    p5, p6, p7, p8 = p6, p5, p8, p7
                elif p5 == p7:
                    pass
                else:
                    continue
                
                if p2 != p6:
                    p5, p6, p7, p8 = p7, p8, p5, p6

                if p1 != p5 or p2 != p6 or p4 == p8:
                    continue

                col = next((c for c in collinears if set({p1, p2, p4}) == set(c.points)), None)
                if not col:
                    p1, p2, p3, p4, p5, p6, p7, p8 = p5, p6, p7, p8, p1, p2, p3, p4
                    col = next((c for c in collinears if set({p1, p2, p4}) == set(c.points)), None)
                    if not col:
                        continue

                new_relations.append(Perpendicular(
                    p2, p8, p4, p8,
                    parents=[cong1, cong2, col],
                    rule="congOAOB_congOBOC_colOAB__perpACBC"
                ))

        return new_relations

    def congOAOB_congOAOC_congOAOD__cyclicABCD(self) -> List[RelationNode]:
        new_relations = []
        congruences = self.problem.relations.get("cong", [])
        for i in range(len(congruences)):
            cong1 = congruences[i]
            p1, p2, p3, p4 = cong1.points
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

            seg1 = frozenset({p1, p2})
            congruent_segs = []
            for j in range(len(congruences)):
                if i == j:
                    continue
                p5, p6, p7, p8 = congruences[j].points
                if frozenset({p5, p6}) == seg1:
                    if p1 == p8:
                        p7, p8 = p8, p7
                    elif p1 == p7:
                        pass
                    else:
                        continue
                    congruent_segs.append(((p7, p8), congruences[j]))
                elif frozenset({p7, p8}) == seg1:
                    if p1 == p6:
                        p5, p6 = p6, p5
                    elif p1 == p5:
                        pass
                    else:
                        continue
                    congruent_segs.append(((p5, p6), congruences[j]))
                
            for j in range(len(congruent_segs)):
                for k in range(j + 1, len(congruent_segs)):
                    (p5, p6), cong2 = congruent_segs[j]
                    (p7, p8), cong3 = congruent_segs[k]
                    if len({p2, p4, p6, p8}) != 4:
                        continue

                    new_relations.append(Cyclic(
                        p2, p4, p6, p8,
                        parents=[cong1, cong2, cong3],
                        rule="congOAOB_congOAOC_congOAOD__cyclicABCD"
                    ))

            p1, p2, p3, p4 = p3, p4, p1, p2
            seg1 = frozenset({p1, p2})
            congruent_segs = []
            for j in range(len(congruences)):
                if i == j:
                    continue
                p5, p6, p7, p8 = congruences[j].points
                if frozenset({p5, p6}) == seg1:
                    if p1 == p8:
                        p7, p8 = p8, p7
                    elif p1 == p7:
                        pass
                    else:
                        continue
                    congruent_segs.append(((p7, p8), congruences[j]))
                elif frozenset({p7, p8}) == seg1:
                    if p1 == p6:
                        p5, p6 = p6, p5
                    elif p1 == p5:
                        pass
                    else:
                        continue
                    congruent_segs.append(((p5, p6), congruences[j]))
                
            for j in range(len(congruent_segs)):
                for k in range(j + 1, len(congruent_segs)):
                    (p5, p6), cong2 = congruent_segs[j]
                    (p7, p8), cong3 = congruent_segs[k]
                    if len({p2, p4, p6, p8}) != 4:
                        continue

                    new_relations.append(Cyclic(
                        p2, p4, p6, p8,
                        parents=[cong1, cong2, cong3],
                        rule="congOAOB_congOAOC_congOAOD__cyclicABCD"
                    ))
            
        return new_relations
                           
    def midpMAB_midpNAC__paraMNBC(self) -> List[RelationNode]:
        new_relations = []
        midpoints = self.problem.relations.get("midp", [])
        for i in range(len(midpoints)):
            for j in range(i + 1, len(midpoints)):
                midp1 = midpoints[i]
                midp2 = midpoints[j]
                p1, p2, p3 = midp1.points
                p4, p5, p6 = midp2.points
                if p2 == p5:
                    pass
                elif p2 == p6:
                    p5, p6 = p6, p5
                elif p3 == p5:
                    p2, p3 = p3, p2
                elif p3 == p6:
                    p2, p3, p5, p6 = p3, p2, p6, p5
                else:
                    continue

                if p1 == p4 or p3 == p6:
                    continue

                new_relations.append(Parallel(
                    p1, p4, p3, p6,
                    parents=[midp1, midp2],
                    rule="midpMAB_midpNAC__paraMNBC"
                ))

        return new_relations
    
    """
    Helper methods for specific rules starts here.
    """
    
    # checks sameclock
    def is_sameclock(self, p1: Point, p2: Point, p3: Point, p4: Point, p5: Point, p6: Point) -> bool:
        return ((p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x)) * \
                ((p5.x - p4.x) * (p6.y - p4.y) - (p5.y - p4.y) * (p6.x - p4.x)) > 0
    
    # checks if cong relation exists, returns (bool, Congruent or None)
    def is_cong(self, p1: Point, p2: Point, p3: Point, p4: Point) -> Tuple[bool, Optional[Congruent]]:
        if set({p1, p2}) == set({p3, p4}):
            return True, None
        else:
            cong_rels = self.problem.relations.get("cong", [])
            for cong in cong_rels:
                if cong.relation == frozenset({frozenset({p1, p2}), frozenset({p3, p4})}):
                    return True, cong
                
        return False, None
