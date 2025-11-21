from matplotlib import pyplot as plt
from typing import Dict, Tuple
import numpy as np
from relations import *

plt.style.use("seaborn-v0_8-whitegrid")  # clean white background


class Canva:
    def __init__(self, points: List[RelationNode],
                 points_dict: Dict[str, Tuple[float, float]], 
                 lines: Dict[str, Tuple[float, float, float]], 
                 circles: Dict[str, Tuple[float, float, float]]):
        self.fig, self.ax = setup_geometry_plot()
        self.points = points
        self.points_dict = points_dict
        self.auxiliary_points = []
        self.auxiliary_points_dict = {}
        self.lines = lines
        self.circles = circles
        self.auxiliary_counter = 1

    def plot(self):
        plot_points(self.ax, self.points_dict)
        plot_points(self.ax, self.auxiliary_points_dict)
        plot_lines_from_eq(self.ax, self.lines)
        plot_circles_from_eq(self.ax, self.circles)
        plt.show()

    def add_point(self, x: float, y: float) -> Point:
        name = "X" + str(self.auxiliary_counter)
        self.auxiliary_counter += 1
        p = Point(name, x, y)
        self.auxiliary_points.append(p)
        self.auxiliary_points_dict[name] = (x, y)
        return p

    def free(self) -> Point:
        x = np.random.uniform(-5, 5)
        y = np.random.uniform(-5, 5)
        return self.add_point(x, y)
    
    def midpoint(self, a: Point, b: Point) -> Tuple[Point, List[RelationNode]]:
        x = (a.x + b.x) / 2
        y = (a.y + b.y) / 2
        p = self.add_point(x, y)
        midp = Midpoint(p, a, b, rule="construction")
        return p, [midp]
    
    def mirror(self, a: Point, b: Point) -> Tuple[Point, List[RelationNode]]:
        # mirror point a across point b
        x = 2 * b.x - a.x
        y = 2 * b.y - a.y
        p = self.add_point(x, y)
        midp = Midpoint(b, a, p, rule="construction")
        return p, [midp]

    def intersect_lines(self, a: Point, b: Point, c: Point, d: Point) -> Tuple[Point, List[RelationNode]]:
        if isinstance(a, np.ndarray):
            a = Point("temp_a", a[0], a[1])
        if isinstance(b, np.ndarray):
            b = Point("temp_b", b[0], b[1])
        if isinstance(c, np.ndarray):
            c = Point("temp_c", c[0], c[1])
        if isinstance(d, np.ndarray):
            d = Point("temp_d", d[0], d[1])
        if is_collinear(a, b, c) or is_collinear(a, b, d) or is_collinear(c, d, a) or is_collinear(c, d, b):
            return None, []
        m1 = (a.y - b.y) / (a.x - b.x)
        m2 = (c.y - d.y) / (c.x - d.x)
        c1 = a.y - m1 * a.x
        c2 = c.y - m2 * c.x
        x = (c1 - c2) / (m2 - m1)
        y = x * m1 + c1
        p = self.add_point(x, y)
        col1 = Collinear(a, b, p, rule="construction")
        col2 = Collinear(c, d, p, rule="construction")
        return p, [col1, col2]

    def anglebisector(self, a: Point, b: Point, c: Point) -> Tuple[Point, List[RelationNode]]:
        # bisects angle ABC, intersects on AC at X
        if is_collinear(a, b, c):
            return None, []
        A, B, C = np.array([a.x, a.y]), np.array([b.x, b.y]), np.array([c.x, c.y])
        BA, BC = A - B, C - B
        BA_norm = BA / np.linalg.norm(BA)
        BC_norm = BC / np.linalg.norm(BC)
        bisector_dir = (BA_norm + BC_norm) / np.linalg.norm(BA_norm + BC_norm)
        start_point = B - bisector_dir * 100
        end_point = B + bisector_dir * 100
        X, _ = self.intersect_lines(start_point, end_point, a, c)
        if X is None:
            return None, []
        eqangle = EqualAngle(a, b, X, X, b, c, rule="construction")
        col = Collinear(a, c, X, rule="construction")
        line_name = f"Line_{b.name}{X.name}"
        self.lines[line_name] = (bisector_dir[1], -bisector_dir[0], bisector_dir[0] * B[1] - bisector_dir[1] * B[0])
        return X, [eqangle, col]

    def foot(self, pt_from: Point, a: Point, b: Point) -> Tuple[Point, List[RelationNode]]:
        if is_collinear(pt_from, a, b):
            return None, []
        A, B = np.array([a.x, a.y]), np.array([b.x, b.y])
        AB = A - B
        perp = np.array([-AB[1], AB[0]])
        perp /= np.linalg.norm(perp)
        start_point = np.array([pt_from.x, pt_from.y]) - perp * 100
        end_point = np.array([pt_from.x, pt_from.y]) + perp * 100
        X, _ = self.intersect_lines(start_point, end_point, a, b)
        if X is None:
            return None, []
        perpendicular = Perpendicular(pt_from, X, a, b, rule="construction")
        col = Collinear(X, a, b, rule="construction")
        line_name = f"Line_{pt_from.name}{X.name}"
        # use start and end points to define line equation
        self.lines[line_name] = (perp[1], -perp[0], perp[0] * start_point[1] - perp[1] * start_point[0])
        return X, [perpendicular, col]

    def circle(self, a: Point, b: Point, c: Point) -> Tuple[Point, List[RelationNode]]:
        # find using perp bisectors of 2 sides intersecting
        if is_collinear(a, b, c):
            return None, []
        midp_AB_x = (a.x + b.x) / 2
        midp_AB_y = (a.y + b.y) / 2
        midp_BC_x = (b.x + c.x) / 2
        midp_BC_y = (b.y + c.y) / 2
        A, B = np.array([a.x, a.y]), np.array([b.x, b.y])
        AB = A - B
        perp = np.array([-AB[1], AB[0]])
        perp /= np.linalg.norm(perp)
        AB_start_point = np.array([midp_AB_x, midp_AB_y]) - perp * 100
        AB_end_point = np.array([midp_AB_x, midp_AB_y]) + perp * 100
        B, C = np.array([b.x, b.y]), np.array([c.x, c.y])
        BC = B - C
        perp = np.array([-BC[1], BC[0]])
        perp /= np.linalg.norm(perp)
        BC_start_point = np.array([midp_BC_x, midp_BC_y]) - perp * 100
        BC_end_point = np.array([midp_BC_x, midp_BC_y]) + perp * 100
        X, _ = self.intersect_lines(AB_start_point, AB_end_point, BC_start_point, BC_end_point)
        if X is None:
            return None, []
        circle = Circle(X, a, b, c, rule="construction")
        circle_name = f"Circle_{X.name}{a.name}{b.name}{c.name}"
        self.circles[circle_name] = (X.x, X.y, np.sqrt((a.x - X.x)**2 + (a.y - X.y)**2))
        return X, [circle]
    
    def incenter(self, a: Point, b: Point, c: Point) -> Tuple[List[Point], List[RelationNode]]:
        # find the intersection of 2 angle bisectors
        if is_collinear(a, b, c):
            return [], []
        A, B, C = np.array([a.x, a.y]), np.array([b.x, b.y]), np.array([c.x, c.y])
        BA, BC = A - B, C - B
        BA_norm = BA / np.linalg.norm(BA)
        BC_norm = BC / np.linalg.norm(BC)
        bisector_dir = (BA_norm + BC_norm) / np.linalg.norm(BA_norm + BC_norm)
        B_start_point = B - bisector_dir * 100
        B_end_point = B + bisector_dir * 100
        AB = B - A
        AC = C - A
        AB_norm = AB / np.linalg.norm(AB)
        AC_norm = AC / np.linalg.norm(AC)
        bisector_dir = (AB_norm + AC_norm) / np.linalg.norm(AB_norm + AC_norm)
        A_start_point = A - bisector_dir * 100
        A_end_point = A + bisector_dir * 100
        I, _ = self.intersect_lines(A_start_point, A_end_point, B_start_point, B_end_point)
        X, rel1 = self.anglebisector(b, a, c)
        Y, rel2 = self.anglebisector(c, b, a)
        Z, rel3 = self.anglebisector(a, c, b)
        relations = rel1 + rel2 + rel3
        if I is None or X is None or Y is None or Z is None:
            return [], []
        eqangle1 = EqualAngle(a, b, I, I, b, c, rule="construction")
        eqangle2 = EqualAngle(b, c, I, I, c, a, rule="construction")
        eqangle3 = EqualAngle(c, a, I, I, a, b, rule="construction")
        col1 = Collinear(a, I, X, rule="construction")
        col2 = Collinear(b, I, Y, rule="construction")
        col3 = Collinear(c, I, Z, rule="construction")
        line1_name = f"Line_{a.name}{I.name}{X.name}"
        line2_name = f"Line_{b.name}{I.name}{Y.name}"
        line3_name = f"Line_{c.name}{I.name}{Z.name}"
        self.lines[line1_name] = (bisector_dir[1], -bisector_dir[0], bisector_dir[0] * A[1] - bisector_dir[1] * A[0])
        self.lines[line2_name] = (bisector_dir[1], -bisector_dir[0], bisector_dir[0] * B[1] - bisector_dir[1] * B[0])
        self.lines[line3_name] = (bisector_dir[1], -bisector_dir[0], bisector_dir[0] * C[1] - bisector_dir[1] * C[0])
        return [I, X, Y, Z], relations + [eqangle1, eqangle2, eqangle3, col1, col2, col3]
    
    def incenter2(self, a: Point, b: Point, c: Point) -> Tuple[List[Point], List[RelationNode]]:
        points, relations = self.incenter(a, b, c)
        if not points:
            return [], []
        I, X, Y, Z = points
        X1, rel1 = self.foot(I, b, c)
        Y1, rel2 = self.foot(I, a, c)
        Z1, rel3 = self.foot(I, a, b)
        relations += rel1 + rel2 + rel3
        if X1 is None or Y1 is None or Z1 is None:
            return [], []
        perp1 = Perpendicular(a, b, I, Z1, rule="construction")
        perp2 = Perpendicular(a, c, I, Y1, rule="construction")
        perp3 = Perpendicular(b, c, I, X1, rule="construction")
        circle = Circle(I, X1, Y1, Z1, rule="construction")
        cong = Congruent(I, X1, I, X, rule="construction")
        circle_name = f"Circle_{I.name}{X1.name}{Y1.name}{Z1.name}"
        self.circles[circle_name] = (I.x, I.y, np.sqrt((X1.x - I.x)**2 + (X1.y - I.y)**2))
        return [I, X, Y, Z, X1, Y1, Z1], relations + [perp1, perp2, perp3, circle, cong]
    
    def excenter(self, a: Point, b: Point, c: Point) -> Tuple[List[Point], List[RelationNode]]:
        # find the intersection of internal bisector of angle A and external bisectors of angles B and C
        if is_collinear(a, b, c):
            return [], []
        A, B, C = np.array([a.x, a.y]), np.array([b.x, b.y]), np.array([c.x, c.y])
        AB = B - A
        AC = C - A
        AB_norm = AB / np.linalg.norm(AB)
        AC_norm = AC / np.linalg.norm(AC)
        bisector_dir_A = (AB_norm + AC_norm) / np.linalg.norm(AB_norm + AC_norm)
        A_start_point = A - bisector_dir_A * 100
        A_end_point = A + bisector_dir_A * 100
        BA, BC = A - B, C - B
        BA_norm = BA / np.linalg.norm(BA)
        BC_norm = BC / np.linalg.norm(BC)
        bisector_dir_B = (BA_norm - BC_norm) / np.linalg.norm(BA_norm - BC_norm)
        B_start_point = B - bisector_dir_B * 100
        B_end_point = B + bisector_dir_B * 100
        CA = A - C
        CB = B - C
        CA_norm = CA / np.linalg.norm(CA)
        CB_norm = CB / np.linalg.norm(CB)
        bisector_dir_C = (CA_norm - CB_norm) / np.linalg.norm(CA_norm - CB_norm)
        C_start_point = C - bisector_dir_C * 100
        C_end_point = C + bisector_dir_C * 100

        I, _ = self.intersect_lines(A_start_point, A_end_point, B_start_point, B_end_point)
        if I is None:
            return [], []
        X, _ = self.intersect_lines(A_start_point, A_end_point, b, c)  # on BC
        Y, _ = self.intersect_lines(B_start_point, B_end_point, a, c)  # on AC
        Z, _ = self.intersect_lines(C_start_point, C_end_point, a, b)  # on AB
        if X is None or Y is None or Z is None:
            return [], []
        
        eqangle1 = EqualAngle(b, a, I, I, a, c, rule="construction")
        eqangle2 = EqualAngle(c, b, I, I, b, a, rule="construction")
        eqangle3 = EqualAngle(a, c, I, I, c, b, rule="construction")
        col1 = Collinear(a, X, I, rule="construction")
        col2 = Collinear(b, Y, I, rule="construction")
        col3 = Collinear(c, Z, I, rule="construction")
        col4 = Collinear(b, c, X, rule="construction")
        col5 = Collinear(a, c, Y, rule="construction")
        col6 = Collinear(a, b, Z, rule="construction")

        line1_name = f"Line_{a.name}{I.name}{X.name}"
        self.lines[line1_name] = (bisector_dir_A[1], -bisector_dir_A[0], bisector_dir_A[0] * A[1] - bisector_dir_A[1] * A[0])
        line2_name = f"Line_{b.name}{I.name}{Y.name}"
        self.lines[line2_name] = (bisector_dir_B[1], -bisector_dir_B[0], bisector_dir_B[0] * B[1] - bisector_dir_B[1] * B[0])
        line3_name = f"Line_{c.name}{I.name}{Z.name}"
        self.lines[line3_name] = (bisector_dir_C[1], -bisector_dir_C[0], bisector_dir_C[0] * C[1] - bisector_dir_C[1] * C[0])
        
        relations = [eqangle1, eqangle2, eqangle3, col1, col2, col3, col4, col5, col6]
        return [I, X, Y, Z], relations
    
    def excenter2(self, a: Point, b: Point, c: Point) -> Tuple[List[Point], List[RelationNode]]:
        points, relations = self.excenter(a, b, c)
        if not points:
            return [], []
        I, X, Y, Z = points
        X1, rel1 = self.foot(I, b, c)
        Y1, rel2 = self.foot(I, a, c)
        Z1, rel3 = self.foot(I, a, b)
        relations += rel1 + rel2 + rel3
        if X1 is None or Y1 is None or Z1 is None:
            return [], []
        perp1 = Perpendicular(a, b, I, Z1, rule="construction")
        perp2 = Perpendicular(a, c, I, Y1, rule="construction")
        perp3 = Perpendicular(b, c, I, X1, rule="construction")
        circle = Circle(I, X1, Y1, Z1, rule="construction")
        cong = Congruent(I, X1, I, X, rule="construction")
        circle_name = f"Circle_{I.name}{X1.name}{Y1.name}{Z1.name}"
        self.circles[circle_name] = (I.x, I.y, np.sqrt((X1.x - I.x)**2 + (X1.y - I.y)**2))
        return [I, X, Y, Z, X1, Y1, Z1], relations + [perp1, perp2, perp3, circle, cong]

    def centroid(self, a: Point, b: Point, c: Point) -> Tuple[List[Point], List[RelationNode]]:
        if is_collinear(a, b, c):
            return [], []
        x = (a.x + b.x + c.x) / 3
        y = (a.y + b.y + c.y) / 3
        p = self.add_point(x, y)
        X, midp1 = self.midpoint(b, c)
        Y, midp2 = self.midpoint(a, c)
        Z, midp3 = self.midpoint(a, b)
        col1 = Collinear(a, X, p, rule="construction")
        col2 = Collinear(b, Y, p, rule="construction")
        col3 = Collinear(c, Z, p, rule="construction")
        line1_name = f"Line_{a.name}{X.name}"
        line2_name = f"Line_{b.name}{Y.name}"
        line3_name = f"Line_{c.name}{Z.name}"
        self.lines[line1_name] = (X.y - a.y, a.x - X.x, X.x * a.y - X.y * a.x)
        self.lines[line2_name] = (Y.y - b.y, b.x - Y.x, Y.x * b.y - Y.y * b.x)
        self.lines[line3_name] = (Z.y - c.y, c.x - Z.x, Z.x * c.y - Z.y * c.x)
        return [p, X, Y, Z], midp1 + midp2 + midp3 + [col1, col2, col3]
    
    def orthocenter(self, a: Point, b: Point, c: Point) -> Tuple[Point, List[RelationNode]]:
        if is_collinear(a, b, c):
            return None, []
        A, B, C = np.array([a.x, a.y]), np.array([b.x, b.y]), np.array([c.x, c.y])
        AB = B - A
        AC = C - A
        perp_A = np.array([-AB[1], AB[0]])
        perp_A /= np.linalg.norm(perp_A)
        A_start_point = A - perp_A * 100
        A_end_point = A + perp_A * 100
        BC = C - B
        perp_B = np.array([-BC[1], BC[0]])
        perp_B /= np.linalg.norm(perp_B)
        B_start_point = B - perp_B * 100
        B_end_point = B + perp_B * 100
        H, _ = self.intersect_lines(A_start_point, A_end_point, B_start_point, B_end_point)
        if H is None:
            return None, []
        perp1 = Perpendicular(a, b, H, c, rule="construction")
        perp2 = Perpendicular(b, c, H, a, rule="construction")
        perp3 = Perpendicular(c, a, H, b, rule="construction")
        return H, [perp1, perp2, perp3]
    
    def orthocenter2(self, a: Point, b: Point, c: Point) -> Tuple[List[Point], List[RelationNode]]:
        H, relations = self.orthocenter(a, b, c)
        if H is None:
            return [], []
        X, rel1 = self.foot(a, b, c)
        Y, rel2 = self.foot(b, a, c)
        Z, rel3 = self.foot(c, a, b)
        relations += rel1 + rel2 + rel3
        if X is None or Y is None or Z is None:
            return [], []
        col1 = Collinear(a, H, X, rule="construction")
        col2 = Collinear(b, H, Y, rule="construction")
        col3 = Collinear(c, H, Z, rule="construction")
        return [H, X, Y, Z], relations + [col1, col2, col3]
    
    def reflect(self, a: Point, b: Point, c: Point) -> Tuple[List[Point], List[RelationNode]]:
        # reflect point a across line BC
        Y, rel1 = self.foot(a, b, c)
        if Y is None:
            return None, []
        X, rel2 = self.mirror(a, Y)
        if X is None:
            return None, []
        col = Collinear(Y, b, c, rule="construction")
        perp = Perpendicular(a, X, b, c, rule="construction")
        return [X, Y], rel1 + rel2 + [col, perp]

    def parallel(self, a: Point, b: Point, c: Point) -> Tuple[Point, List[RelationNode]]:
        # line through point a parallel to line bc
        A, B, C = np.array([a.x, a.y]), np.array([b.x, b.y]), np.array([c.x, c.y])
        BC = C - B
        BC_norm = BC / np.linalg.norm(BC)
        parallel_dir = BC_norm
        line_name = f"Line_{a.name}_parallel_{b.name}{c.name}"
        self.lines[line_name] = (parallel_dir[1], -parallel_dir[0], parallel_dir[0] * A[1] - parallel_dir[1] * A[0])
        D = A + parallel_dir * np.linalg.norm(BC)
        d = self.add_point(D[0], D[1])
        para = Parallel(a, d, b, c, rule="construction")
        cong = Congruent(b, c, a, d, rule="construction")
        return d, [para, cong]

    def tangent(self, a: Point, b: Point) -> Tuple[Point, List[RelationNode]]:
        A, B = np.array([a.x, a.y]), np.array([b.x, b.y])
        AB = B - A
        perp = np.array([-AB[1], AB[0]])
        perp /= np.linalg.norm(perp)
        C = B + perp * np.linalg.norm(AB)
        c = self.add_point(C[0], C[1])
        line_name = f"Line_{b.name}{c.name}"
        self.lines[line_name] = (perp[1], -perp[0], perp[0] * B[1] - perp[1] * B[0])
        perp = Perpendicular(a, b, b, c, rule="construction")
        cong = Congruent(a, b, b, c, rule="construction")
        return c, [perp, cong]

    def tangent2(self, a: Point, o: Point, b: Point) -> Tuple[List[Point], List[RelationNode]]:
        # tangent from point a to circle with center o passing through b
        A, O, B = np.array([a.x, a.y]), np.array([o.x, o.y]), np.array([b.x, b.y])
        OB = B - O
        r = np.linalg.norm(OB)
        AO = O - A
        d = np.linalg.norm(AO)
        if d - r < 1e-10:
            return [], []
        l = np.sqrt(d**2 - r**2)
        AO_norm = AO / d
        perp_dir = np.array([-AO_norm[1], AO_norm[0]])
        T1 = A + AO_norm * (r**2 / d) + perp_dir * (r * l / d)
        T2 = A + AO_norm * (r**2 / d) - perp_dir * (r * l / d)
        t1 = self.add_point(T1[0], T1[1])
        t2 = self.add_point(T2[0], T2[1])
        line_name1 = f"Line_{a.name}{t1.name}"
        line_name2 = f"Line_{a.name}{t2.name}"
        circle_name = f"Circle_{o.name}{b.name}{t1.name}{t2.name}"
        self.lines[line_name1] = (t1.y - a.y, a.x - t1.x, t1.x * a.y - t1.y * a.x)
        self.lines[line_name2] = (t2.y - a.y, a.x - t2.x, t2.x * a.y - t2.y * a.x)
        self.circles[circle_name] = (o.x, o.y, r)
        circle = Circle(o, b, t1, t2, rule="construction")
        perp1 = Perpendicular(a, t1, o, t1, rule="construction")
        perp2 = Perpendicular(a, t2, o, t2, rule="construction")
        return [t1, t2], [circle, perp1, perp2]
    
    def on_dia(self, a: Point, b: Point) -> Tuple[Point, List[RelationNode]]:
        # create a point c such that perp a c c b
        A, B = np.array([a.x, a.y]), np.array([b.x, b.y])
        O = np.array([(a.x + b.x) / 2, (a.y + b.y) / 2])
        r = np.linalg.norm(A - O)
        theta = np.random.uniform(0, 2 * np.pi)
        Cx = O[0] + r * np.cos(theta)
        Cy = O[1] + r * np.sin(theta)
        c = self.add_point(Cx, Cy)
        perp = Perpendicular(a, c, c, b, rule="construction")
        line_name1 = f"Line_{a.name}{c.name}"
        line_name2 = f"Line_{c.name}{b.name}"
        self.lines[line_name1] = (c.y - a.y, a.x - c.x, c.x * a.y - c.y * a.x)
        self.lines[line_name2] = (b.y - c.y, c.x - b.x, b.x * c.y - b.y * c.x)
        return c, [perp]


def setup_geometry_plot(square: bool = True):
    fig, ax = plt.subplots(figsize=(8, 8))
    if square:
        ax.set_aspect('equal', 'box')

    # Square gridlines
    ax.grid(True, which='both', linestyle='--', color='lightgray', alpha=0.7)
    ax.minorticks_on()
    ax.tick_params(which='major', length=5, width=0.8, color='dimgray')
    ax.tick_params(which='minor', length=2.5, width=0.5, color='lightgray')

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('lightgray')
        spine.set_linewidth(1)

    ax.set_facecolor('white')
    ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax.yaxis.set_major_locator(plt.MultipleLocator(1))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))

    return fig, ax


def plot_points(ax, points: Dict[str, Tuple[float, float]]):
    for name, (x, y) in points.items():
        ax.plot(x, y, 'o', color='black', markersize=7, zorder=5)
        ax.text(x + 0.05, y + 0.05, name, fontsize=10, color='black', weight='bold')


def plot_lines_from_eq(ax, lines: Dict[str, Tuple[float, float, float]]):
    for line, (a, b, c) in lines.items():
        x_vals = np.linspace(*ax.get_xlim(), 300)
        y_vals = (-a * x_vals - c) / b
        ax.plot(x_vals, y_vals, color='black', lw=1.0, alpha=0.9)


def plot_circles_from_eq(ax, circles: Dict[str, Tuple[float, float, float]]):
    for name, (x, y, r) in circles.items():
        theta = np.linspace(0, 2 * np.pi, 400)
        x_coords = x + r * np.cos(theta)
        y_coords = y + r * np.sin(theta)
        ax.plot(x_coords, y_coords, color='black', lw=1.0, alpha=0.9)
        ax.text(x, y, "", fontsize=9, color='black', ha='center', va='center')


def is_collinear(a: Point, b: Point, c: Point, tol: float = 1e-10) -> bool:
    area = a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)
    return abs(area) < tol
