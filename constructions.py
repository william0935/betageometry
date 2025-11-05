from matplotlib import pyplot as plt
from typing import Dict, Tuple
import numpy as np
from relations import *
from ar import *
from dd_ar import *

plt.style.use("seaborn-v0_8-whitegrid")  # clean white background

class ConstructionPlotter:
    counter = 1
    plotted_auxillaries = {} # these plotted auxillary points will be of the form X_N for some integer N


    def plot_points(self, ax, points: Dict[str, Tuple[float, float]]):
        for name, (x, y) in points.items():
            ax.plot(x, y, 'o', color='black', markersize=7, zorder=5)
            ax.text(x + 0.05, y + 0.05, name, fontsize=10, color='black', weight='bold')

    def plot_lines_from_eq(self, ax, lines: Dict[str, Tuple[float, float, float]]):
        for line, (a, b, c) in lines.items():
            x_vals = np.linspace(*ax.get_xlim(), 300)
            y_vals = (-a * x_vals - c) / b
            ax.plot(x_vals, y_vals, color='black', lw=1.0, alpha=0.9)

    def plot_circles_from_eq(self, ax, circles: Dict[str, Tuple[float, float, float]]):
        for name, (x, y, r) in circles.items():
            theta = np.linspace(0, 2 * np.pi, 400)
            x_coords = x + r * np.cos(theta)
            y_coords = y + r * np.sin(theta)
            ax.plot(x_coords, y_coords, color='black', lw=1.0, alpha=0.9)
            ax.text(x, y, name, fontsize=9, color='black', ha='center', va='center')

    def plot_newpoint(self, ax, point: Point):
        ax.plot(point.x, point.y, 'o', color='black', markersize=7, zorder=5)
        ax.text(point.x + 0.05, point.y + 0.05, point.name, fontsize=10, color='black', weight='bold')

    def plot_newline(self, ax, a: Point, b: Point):
        ax.plot([a.x, b.x], [a.y, b.y], color='black', lw=1.0, alpha=0.9, linestyle='-')

    def plot_newcircle(self, ax, center: Point, point: Point):
        theta = np.linspace(0, 2 * np.pi, 400)
        r = math.dist(np.array([center.x, center.y]), np.array([point.x, point.y]))
        x_coords = center.x + r * np.cos(theta)
        y_coords = center.y + r * np.sin(theta)
        ax.plot(x_coords, y_coords, color='black', lw=1.0, alpha=0.9)
        ax.text(center.x, center.y, center.name, fontsize=9, color='black', ha='center', va='center')
        
    def plot_midpoint(self, ax, name: str, a: Point, b: Point):
        new_x = (a.x + b.x) / 2
        new_y = (a.y + b.y) / 2
        X = Point(new_x, new_y)
        ax.plot(new_x, new_y, 'o', color='black', markersize=6, zorder=5)
        ax.text(new_x + 0.05, new_y + 0.05, name, fontsize=9, color='black', weight='bold')
        return X



    # bisects angle ABC, intersects on AC at X
    def plot_anglebisector(self, ax, dd_obj: DDWithAR, a: Point, b: Point, c: Point):
        A, B, C = np.array([a.x, a.y]), np.array([b.x, b.y]), np.array([c.x, c.y])
        BA, BC = A - B, C - B
        BA_norm = BA / np.linalg.norm(BA)
        BC_norm = BC / np.linalg.norm(BC)
        bisector_dir = (BA_norm + BC_norm) / np.linalg.norm(BA_norm + BC_norm)
        start_point = B - bisector_dir * 100
        end_point = B + bisector_dir * 100
        X = self.find_intersection(start_point, end_point, a, c)
        ax.plot([start_point[0], end_point[0]], [start_point[1], end_point[1]],
                color='black', lw=1.0, linestyle='--', alpha=0.9)
        ax.plot(X.x, X.y, 'o', color='black', markersize=7, zorder=5)
        ax.text(X.x + 0.05, X.y + 0.05, "X", fontsize=10, color='black', weight='bold')
        new_relations = [EqualAngle(a, b, X, X, b, c), Collinear(a, c, X)]
        a_name = a.name
        b_name = b.name
        c_name = c.name
        dd_obj.angle_table.add_col(a_name + "X")
        dd_obj.angle_table.add_col(b_name + "X")
        dd_obj.angle_table.add_col(c_name + "X")
        for relation in new_relations:
            dd_obj.update_AR_tables_with_relation(relation)
        return new_relations


    def plot_perp(self, ax, dd_obj: DDWithAR, a: Point, b: Point, pt_from: Point):
        A, B = np.array([a.x, a.y]), np.array([b.x, b.y])
        AB = A - B
        perp = np.array([-AB[1], AB[0]])
        perp /= np.linalg.norm(perp)
        start_point = np.array([pt_from.x, pt_from.y]) - perp * 100
        end_point = np.array([pt_from.x, pt_from.y]) + perp * 100
        X = self.find_intersection(a, b, start_point, end_point)
        new_relations = [Perpendicular(pt_from, X, a, b), Collinear(X, a, b)]
        a_name = a.name
        b_name = b.name
        pt_from_name = pt_from.name
        dd_obj.angle_table.add_col(a.name + "X")
        dd_obj.angle_table.add_col(b.name + "X")
        for relation in new_relations:
            dd_obj.update_AR_tables_with_relation(relation)
        ax.plot([start_point[0], end_point[0]], [start_point[1], end_point[1]],
                color='black', lw=1.0, linestyle=':', alpha=0.9)
        return new_relations
        
    def circumcenter(self, ax, dd_obj: DDWithAR, A: Point, B: Point, C: Point):
        # find using perp bisectors of 2 sides intersecting
        midp_AB_x = (A.x + B.x) / 2
        midp_AB_y = (A.y + B.y) / 2
        midp_AB = Point(midp_AB_x, midp_AB_y)
        midp_BC_x = (B.x + C.x) / 2
        midp_BC_y = (B.y + C.y) / 2
        midp_BC = Point(midp_BC_x, midp_BC_y)
        a, b = np.array([A.x, A.y]), np.array([B.x, B.y])
        ab = a - b
        perp = np.array([-ab[1], ab[0]])
        perp /= np.linalg.norm(perp)
        AB_start_point = np.array([midp_AB.x, midp_AB.y]) - perp * 100
        AB_end_point = np.array([midp_AB.x, midp_AB.y]) + perp * 100
        b, c = np.array([B.x, B.y]), np.array([C.x, C.y])
        bc = b - c
        perp = np.array([-bc[1], bc[0]])
        perp /= np.llinalg.norm(perp)
        BC_start_point = np.array([midp_BC.x, midp_BC.y]) - perp * 100
        BC_end_point = np.array([midp_BC.x, midp_BC.y]) + perp * 100
        X = self.find_intersection(AB_start_point, AB_end_point, BC_start_point, BC_end_point)
        new_relations = [Circle(X, A, B, C), Congruent(X, A, X, B), Congruent(X, A, X, C)]
        A_name = A.name
        B_name = B.name
        C_name = C.name
        dd_obj.angle_table.add_col(A_name + "X")
        dd_obj.angle_table.add_col(B_name + "X")
        dd_obj.angle_table.add_col(C_name + "X")
        for relation in new_relations:
            dd_obj.update_AR_tables_with_relation(relation)
        self.plot_newcircle(ax, X, A)
        return new_relations

    def incenter(self, ax, dd_obj: DDWithAR, A: Point, B: Point, C: Point):
        # find the intersection of 2 angle bisectors
        a, b, c = np.array([A.x, A.y]), np.array([B.x, B.y]), np.array([C.x, C.y])
        ba, bc = a - b, c - b
        BA_norm = ba / np.linalg.norm(ba)
        BC_norm = bc / np.linalg.norm(bc)
        bisector_dir = (BA_norm + BC_norm) / np.linalg.norm(BA_norm + BC_norm)
        B_start_point = B - bisector_dir * 100
        B_end_point = B + bisector_dir * 100
        ac = c - a
        AC_norm = ac / np.linalg.norm(ac)
        bisector_dir = (BA_norm + AC_norm) / np.linalg.norm(BA_norm + AC_norm)
        A_start_point = A - bisector_dir * 100
        A_end_point = A + bisector_dir * 100
        I = self.find_intersection(A_start_point, A_end_point, B_start_point, B_end_point)
        new_relations = [EqualAngle(A, B, I, I, B, C), EqualAngle(B, C, I, I, C, A),
                        EqualAngle(C, A, I, I, A, B)]
        A_name = A.name
        B_name = B.name
        C_name = C.name
        dd_obj.angle_table.add_col(A_name + "I")
        dd_obj.angle_table.add_col(B_name + "I")
        dd_obj.angle_table.add_col(C_name + "I")
        for relation in new_relations:
            dd_obj.update_AR_tables_with_relation(relation)
        self.plot_newpoint(ax, I)
        self.plot_newline(ax, I, A)
        self.plot_newline(ax, I, B)
        self.plot_newline(ax, I, C)
        return new_relations

    def find_intersection(self, a: Point, b: Point, c: Point, d: Point):
        m1 = (a.y-b.y)/(a.x-b.x)
        m2 = (c.y-d.y)/(c.x-d.x)
        c1 = a.y - m1*a.x
        c2 = c.y - m2*c.x
        x_x = (c1-c2)/(m2-m1)
        x_y = x_x*m1 + c1
        X = Point(f"X{self.counter}", x_x, x_y)
        self.counter += 1
        return X


    def setup_geometry_plot(self, square: bool = True):
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


    def equalize_axes(ax, span: float = None):
        """Make x and y axes cover the same total range (e.g., both 10 units)."""
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        x_center = (x_min + x_max) / 2
        y_center = (y_min + y_max) / 2

        if span is None:
            # Automatically choose the larger span to keep all geometry visible
            x_span = x_max - x_min
            y_span = y_max - y_min
            span = max(x_span, y_span)

        half = span / 2
        ax.set_xlim(x_center - half, x_center + half)
        ax.set_ylim(y_center - half, y_center + half)
        ax.set_aspect('equal', 'box')
