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

    def add_point(self, x: float, y: float):
        name = "X" + str(self.auxiliary_counter)
        self.auxiliary_counter += 1
        p = Point(name, x, y)
        self.auxiliary_points.append(p)
        self.auxiliary_points_dict[name] = (x, y)
        return p

    def midpoint(self, a: Point, b: Point) -> Tuple[Point, List[RelationNode]]:
        x = (a.x + b.x) / 2
        y = (a.y + b.y) / 2
        p = self.add_point(x, y)
        midp = Midpoint(p, a, b, rule="construction")
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
        X, relations = self.intersect_lines(start_point, end_point, a, c)
        if X is None:
            return None, []
        eqangle = EqualAngle(a, b, X, X, b, c, rule="construction")
        line_name = f"Line_{b.name}{X.name}"
        self.lines[line_name] = (bisector_dir[1], -bisector_dir[0], bisector_dir[0] * B[1] - bisector_dir[1] * B[0])
        return X, relations + [eqangle]

    def foot(self, pt_from: Point, a: Point, b: Point) -> Tuple[Point, List[RelationNode]]:
        if is_collinear(pt_from, a, b):
            return None, []
        A, B = np.array([a.x, a.y]), np.array([b.x, b.y])
        AB = A - B
        perp = np.array([-AB[1], AB[0]])
        perp /= np.linalg.norm(perp)
        start_point = np.array([pt_from.x, pt_from.y]) - perp * 100
        end_point = np.array([pt_from.x, pt_from.y]) + perp * 100
        X, relations = self.intersect_lines(start_point, end_point, a, b)
        if X is None:
            return None, []
        perpendicular = Perpendicular(pt_from, X, a, b, rule="construction")
        col = Collinear(X, a, b, rule="construction")
        line_name = f"Line_{pt_from.name}{X.name}"
        # use start and end points to define line equation
        self.lines[line_name] = (perp[1], -perp[0], perp[0] * start_point[1] - perp[1] * start_point[0])
        return X, relations + [perpendicular, col]

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
        X, relations = self.intersect_lines(AB_start_point, AB_end_point, BC_start_point, BC_end_point)
        if X is None:
            return None, []
        circle = Circle(X, a, b, c, rule="construction")
        circle_name = f"Circle_{X.name}{a.name}{b.name}{c.name}"
        self.circles[circle_name] = (X.x, X.y, np.sqrt((a.x - X.x)**2 + (a.y - X.y)**2))
        return X, relations + [circle]
    
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
        AC = C - A
        AC_norm = AC / np.linalg.norm(AC)
        bisector_dir = (BA_norm + AC_norm) / np.linalg.norm(BA_norm + AC_norm)
        A_start_point = A - bisector_dir * 100
        A_end_point = A + bisector_dir * 100
        I, relations = self.intersect_lines(A_start_point, A_end_point, B_start_point, B_end_point)
        X, rel1 = self.anglebisector(b, a, c)
        Y, rel2 = self.anglebisector(c, b, a)
        Z, rel3 = self.anglebisector(a, c, b)
        relations += rel1 + rel2 + rel3
        if I is None or X is None or Y is None or Z is None:
            return [], []
        eqangle1 = EqualAngle(a, b, I, I, b, c, rule="construction")
        eqangle2 = EqualAngle(b, c, I, I, c, a, rule="construction")
        eqangle3 = EqualAngle(c, a, I, I, a, b, rule="construction")
        col1 = Collinear(a, b, Z, rule="construction")
        col2 = Collinear(a, c, Y, rule="construction")
        col3 = Collinear(b, c, X, rule="construction")
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
    






    # def circumcenter(ax, x_counter, dd_obj: DDWithAR, A: Point, B: Point, C: Point):
    #     # find using perp bisectors of 2 sides intersecting
    #     midp_AB_x = (A.x + B.x) / 2
    #     midp_AB_y = (A.y + B.y) / 2
    #     midp_AB = Point(midp_AB_x, midp_AB_y)
    #     midp_BC_x = (B.x + C.x) / 2
    #     midp_BC_y = (B.y + C.y) / 2
    #     midp_BC = Point(midp_BC_x, midp_BC_y)
    #     a, b = np.array([A.x, A.y]), np.array([B.x, B.y])
    #     ab = a - b
    #     perp = np.array([-ab[1], ab[0]])
    #     perp /= np.linalg.norm(perp)
    #     AB_start_point = np.array([midp_AB.x, midp_AB.y]) - perp * 100
    #     AB_end_point = np.array([midp_AB.x, midp_AB.y]) + perp * 100
    #     b, c = np.array([B.x, B.y]), np.array([C.x, C.y])
    #     bc = b - c
    #     perp = np.array([-bc[1], bc[0]])
    #     perp /= np.llinalg.norm(perp)
    #     BC_start_point = np.array([midp_BC.x, midp_BC.y]) - perp * 100
    #     BC_end_point = np.array([midp_BC.x, midp_BC.y]) + perp * 100
    #     X = find_intersection(AB_start_point, AB_end_point, BC_start_point, BC_end_point)
    #     X.name = "X" + str(x_counter)
    #     dd_obj.problem.points.append(X)
    #     new_relations = [Circle(X, A, B, C), Congruent(X, A, X, B), Congruent(X, A, X, C)]
    #     dd_obj.angle_table.add_col(A.name + X.name)
    #     dd_obj.angle_table.add_col(B.name + X.name)
    #     dd_obj.angle_table.add_col(C.name + X.name)
    #     dd_obj.ratio_table.add_col(A.name + X.name)
    #     dd_obj.ratio_table.add_col(B.name + X.name)
    #     dd_obj.ratio_table.add_col(C.name + X.name)
    #     for relation in new_relations:
    #         dd_obj.problem.add_relation(relation)
    #         dd_obj.update_AR_tables_with_relation(relation)
    #     plot_newcircle(ax, X, A)
    #     plot_newpoint(ax, X)
    #     return X

    # def excenter(ax, j_counter, dd_obj: DDWithAR, A: Point, B: Point, C: Point):
    #     X = anglebisector(ax, 0, dd_obj, B, A, C)
    #     a, b = np.array([A.x, A.y]), np.array([B.x, B.y])
    #     AB = b - a
    #     AB_norm = AB / np.linalg.norm(AB)
    #     AB_extended = b + AB_norm * 100
    #     extend_AB = Point("ABextension", AB_extended[0], AB_extended[1])
    #     Y = anglebisector(ax, 0, dd_obj, extend_AB, B, C)
    #     J  = find_intersection(A, X, B, Y)
    #     J.name = "J" + str(j_counter)
    #     dd_obj.problem.points.append(J)
    #     new_relations = [EqualAngle(B, A, J, J, A, C)]
    #     dd_obj.angle_table.add_col(A.name + J.name)
    #     for relation in new_relations:
    #         dd_obj.update_AR_tables_with_relation(relation)
    #         dd_obj.problem.add_relation(relation)
    #     plot_newpoint(ax, J)
    #     return J

    # def excenter2(ax, x_counter, y_counter, z_counter, j_counter, dd_obj: DDWithAR, A: Point, B: Point, C: Point):
    #     J = excenter(ax, j_counter, dd_obj, A, B, C)
    #     X = perp(ax, x_counter, dd_obj, B, C, J)
    #     Y = perp(ax, y_counter, dd_obj, A, C, J)
    #     Z = perp(ax, z_counter, dd_obj, A, B, J)
    #     dd_obj.problem.points.extend([X, Y, Z, J])
    #     new_relations = [EqualAngle(B, A, J, J, A, C), EqualAngle(X, B, J, J, B, Z),
    #                     EqualAngle(X, C, J, J, C, Y), EqualAngle(J, X, Y, Z),
    #                     Congruent(J, X, J, Y), Congruent(J, X, J, Z), Perpendicular(B, C, X, J),
    #                     Perpendicular(A, C, Y, J), Perpendicular(A, B, Z, J)]
    #     dd_obj.angle_table.add_col(A.name + J.name)
    #     dd_obj.angle_table.add_col(B.name + J.name)
    #     dd_obj.angle_table.add_col(C.name + J.name)
    #     dd_obj.angle_table.add_col(X.name + B.name)
    #     dd_obj.angle_table.add_col(X.name + C.name)
    #     dd_obj.angle_table.add_col(X.name + J.name)
    #     dd_obj.angle_table.add_col(Y.name + C.name)
    #     dd_obj.angle_table.add_col(Y.name + J.name)
    #     dd_obj.angle_table.add_col(Z.name + B.name)
    #     dd_obj.angle_table.add_col(Z.name + J.name)
    #     dd_obj.ratio_table.add_col(X.name + J.name)
    #     dd_obj.ratio_table.add_col(Y.name + J.name)
    #     dd_obj.ratio_table.add_col(Z.name + J.name)
    #     for relation in new_relations:
    #         dd_obj.problem.add_relation(relation)
    #         dd_obj.update_AR_tables_with_relation(relation)
    #     plot_newpoint(ax, J)
    #     plot_newpoint(ax, X)
    #     plot_newpoint(ax, Y)
    #     plot_newpoint(ax, Z)
    #     plot_newcircle(ax, J, X)
    #     return X, Y, Z, J

    # def centroid(ax, x_counter, y_counter, z_counter, g_counter, dd_obj: DDWithAR, A: Point, B: Point, C: Point):
    #     X = Point("X" + str(x_counter), (B.x + C.x) / 2, (B.y + C.y) / 2)
    #     Y = Point("Y" + str(y_counter), (A.x + C.x) / 2, (A.y + C.y) / 2)
    #     Z = Point("Z" + str(z_counter), (A.x + B.x) / 2, (A.y + B.y) / 2)
    #     G = find_intersection(A, X, B, Y)
    #     G.name = "G" + str(g_counter)
    #     dd_obj.problem.points.extend([X, Y, Z, G])
    #     new_relations = [Collinear(X, B, C), Collinear(Y, A, C), Collinear(Z, A, B),
    #                     Midpoint(X, B, C), Collinear(Y, A, C), Collinear(Z, A, B)]
    #     dd_obj.angle_table.add_col(X.name + B.name)
    #     dd_obj.angle_table.add_col(X.name + C.name)
    #     dd_obj.angle_table.add_col(Y.name + A.name)
    #     dd_obj.angle_table.add_col(Y.name + C.name)
    #     dd_obj.angle_table.add_col(Z.name + A.name)
    #     dd_obj.angle_table.add_col(Z.name + B.name)
    #     dd_obj.ratio_table.add_col(X.name + B.name)
    #     dd_obj.ratio_table.add_col(X.name + C.name)
    #     dd_obj.ratio_table.add_col(Y.name + A.name)
    #     dd_obj.ratio_table.add_col(Y.name + C.name)
    #     dd_obj.ratio_table.add_col(Z.name + A.name)
    #     dd_obj.ratio_table.add_col(Z.name + B.name)
    #     for relation in new_relations:
    #         dd_obj.problem.add_relation(relation)
    #         dd_obj.update_AR_tables_with_relation(relation)
    #     plot_newpoint(ax, G)
    #     plot_newpoint(ax, X)
    #     plot_newpoint(ax, Y)
    #     plot_newpoint(ax, Z)
    #     return G, X, Y, Z













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
