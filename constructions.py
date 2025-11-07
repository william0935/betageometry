from matplotlib import pyplot as plt
from typing import Dict, Tuple
import numpy as np
from relations import *
from ar import *
from dd_ar import *

plt.style.use("seaborn-v0_8-whitegrid")  # clean white background

class ConstructionPlotter:
    counter = 1
    midpoint_counter = 1
    incenter_counter = 1
    orthocenter_counter = 1
    circumcenter_counter = 1
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
        # use tuples for math.dist (works with sequences)
        r = math.dist((center.x, center.y), (point.x, point.y))
        x_coords = center.x + r * np.cos(theta)
        y_coords = center.y + r * np.sin(theta)
        ax.plot(x_coords, y_coords, color='black', lw=1.0, alpha=0.9)

        # Plot the center as a visible point and label it (slight offset)
        ax.plot(center.x, center.y, 'o', color='black', markersize=7, zorder=5)
        ax.text(center.x + 0.05, center.y + 0.05, center.name, fontsize=10, color='black', weight='bold')
        
    def plot_midpoint(self, ax, dd_obj: DDWithAR, a: Point, b: Point):
        new_x = (a.x + b.x) / 2
        new_y = (a.y + b.y) / 2
        X = Point(f"M{self.midpoint_counter}", new_x, new_y)
        ax.plot(new_x, new_y, 'o', color='black', markersize=6, zorder=5)
        ax.text(new_x + 0.05, new_y + 0.05, X.name, fontsize=9, color='black', weight='bold')
        self.plot_newline(ax, a, b)
        self.midpoint_counter += 1
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
        X = self.find_intersection(a, c, start_point, end_point)
        self.counter -= 1
        ax.plot([start_point[0], end_point[0]], [start_point[1], end_point[1]],
                color='black', lw=1.0, linestyle='--', alpha=0.9)
        ax.plot(X.x, X.y, 'o', color='black', markersize=7, zorder=5)
        ax.text(X.x + 0.05, X.y + 0.05, f"X{self.counter}", fontsize=10, color='black', weight='bold')
        self.plot_newline(ax, a, b)
        self.plot_newline(ax, b, c)
        self.plot_newline(ax, a, c)
        new_relations = [EqualAngle(a, b, X, X, b, c), Collinear(a, c, X)]
        a_name = a.name
        b_name = b.name
        c_name = c.name
        dd_obj.angle_table.add_col(a_name + f"X{self.counter}")
        dd_obj.angle_table.add_col(b_name + f"X{self.counter}")
        dd_obj.angle_table.add_col(c_name + f"X{self.counter}")
        for relation in new_relations:
            dd_obj.update_AR_tables_with_relation(relation)
        self.counter += 1
        return X


    def plot_perp(self, ax, dd_obj: DDWithAR, a: Point, b: Point, pt_from: Point):
        A, B = np.array([a.x, a.y]), np.array([b.x, b.y])
        AB = A - B
        perp = np.array([-AB[1], AB[0]])
        perp /= np.linalg.norm(perp)
        start_point = np.array([pt_from.x, pt_from.y]) - perp * 100
        end_point = np.array([pt_from.x, pt_from.y]) + perp * 100
        X = self.find_intersection(a, b, start_point, end_point)
        self.counter -= 1
        new_relations = [Perpendicular(pt_from, X, a, b), Collinear(X, a, b)]
        a_name = a.name
        b_name = b.name
        pt_from_name = pt_from.name
        dd_obj.angle_table.add_col(a.name + f"X{self.counter}")
        dd_obj.angle_table.add_col(b.name + f"X{self.counter}")
        for relation in new_relations:
            dd_obj.update_AR_tables_with_relation(relation)
        ax.plot([start_point[0], end_point[0]], [start_point[1], end_point[1]],
                color='black', lw=1.0, linestyle=':', alpha=0.9)
        self.counter += 1
        return X
        
    def circumcenter(self, ax, dd_obj: DDWithAR, A: Point, B: Point, C: Point):
        """Find circumcenter as intersection of perpendicular bisectors of AB and BC."""
        # points as numpy arrays
        a = np.array([A.x, A.y], dtype=float)
        b = np.array([B.x, B.y], dtype=float)
        c = np.array([C.x, C.y], dtype=float)

        # midpoints as numpy arrays
        mid_AB = (a + b) / 2.0
        mid_BC = (b + c) / 2.0

        # perpendicular direction for AB
        ab = a - b
        perp_ab = np.array([-ab[1], ab[0]], dtype=float)
        perp_ab /= np.linalg.norm(perp_ab)

        # perpendicular direction for BC
        bc = b - c
        perp_bc = np.array([-bc[1], bc[0]], dtype=float)
        perp_bc /= np.linalg.norm(perp_bc)

        # visual ray length (tweak if you want shorter/longer)
        length = 100.0

        AB_start_point = mid_AB - perp_ab * length
        AB_end_point   = mid_AB + perp_ab * length
        BC_start_point = mid_BC - perp_bc * length
        BC_end_point   = mid_BC + perp_bc * length

        # find intersection (find_intersection accepts array-like or Points)
        O = self.find_intersection(AB_start_point, AB_end_point, BC_start_point, BC_end_point)
        self.counter -= 1

        new_relations = [Circle(O, A, B, C), Congruent(O, A, O, B), Congruent(O, A, O, C)]

        A_name = A.name
        B_name = B.name
        C_name = C.name
        dd_obj.angle_table.add_col(A_name + f"O{self.circumcenter_counter}")
        dd_obj.angle_table.add_col(B_name + f"O{self.circumcenter_counter}")
        dd_obj.angle_table.add_col(C_name + f"O{self.circumcenter_counter}")
        for relation in new_relations:
            dd_obj.update_AR_tables_with_relation(relation)

        # draw circle and update counter
        self.plot_newcircle(ax, O, A)
        self.circumcenter_counter += 1
        return O

    def incenter(self, ax, dd_obj: DDWithAR, A: Point, B: Point, C: Point):
        """
        Compute and plot the incenter (intersection of internal angle bisectors).
        Returns a list of the EqualAngle relations (same shape as other construction methods).
        """

        # convert to numpy arrays (vectors)
        a = np.array([A.x, A.y], dtype=float)
        b = np.array([B.x, B.y], dtype=float)
        c = np.array([C.x, C.y], dtype=float)

        # Bisector at A: use vectors AB and AC (both originate at A)
        AB = b - a
        AC = c - a
        ABn = AB / np.linalg.norm(AB)
        ACn = AC / np.linalg.norm(AC)
        bisector_dir_A = (ABn + ACn)
        bisector_dir_A /= np.linalg.norm(bisector_dir_A)

        # Bisector at B: use vectors BA and BC (both originate at B)
        BA = a - b
        BC = c - b
        BAn = BA / np.linalg.norm(BA)
        BCn = BC / np.linalg.norm(BC)
        bisector_dir_B = (BAn + BCn)
        bisector_dir_B /= np.linalg.norm(bisector_dir_B)

        # visual ray length (you can compute from ax.get_xlim()/get_ylim() if you prefer)
        length = 100.0
        A_start = a - bisector_dir_A * length
        A_end   = a + bisector_dir_A * length
        B_start = b - bisector_dir_B * length
        B_end   = b + bisector_dir_B * length

        # find intersection (works with array-like)
        I = self.find_intersection(A_start, A_end, B_start, B_end)

        # build relations (same format you used before)
        new_relations = [
            EqualAngle(A, B, I, I, B, C),
            EqualAngle(B, C, I, I, C, A),
            EqualAngle(C, A, I, I, A, B)
        ]

        # update angle-table columns
        A_name = A.name
        B_name = B.name
        C_name = C.name
        dd_obj.angle_table.add_col(A_name + f"I{self.incenter_counter}")
        dd_obj.angle_table.add_col(B_name + f"I{self.incenter_counter}")
        dd_obj.angle_table.add_col(C_name + f"I{self.incenter_counter}")
        for relation in new_relations:
            dd_obj.update_AR_tables_with_relation(relation)

        # plot the incenter and connecting lines
        self.plot_newpoint(ax, I)
        self.plot_newline(ax, I, A)
        self.plot_newline(ax, I, B)
        self.plot_newline(ax, I, C)

        self.incenter_counter += 1
        return new_relations

    def find_intersection(self, a, b, c, d):
        def _to_xy(p):
            if hasattr(p, "x") and hasattr(p, "y"):
                return float(p.x), float(p.y)
            arr = np.asarray(p, dtype=float)
            return float(arr[0]), float(arr[1])

        x1, y1 = _to_xy(a)
        x2, y2 = _to_xy(b)
        x3, y3 = _to_xy(c)
        x4, y4 = _to_xy(d)

        # Line coefficients for form: A*x + B*y = C
        A1 = y2 - y1
        B1 = x1 - x2
        C1 = A1 * x1 + B1 * y1

        A2 = y4 - y3
        B2 = x3 - x4
        C2 = A2 * x3 + B2 * y3

        det = A1 * B2 - A2 * B1
        if abs(det) < 1e-12:
            raise ValueError("Lines are parallel or nearly parallel; no unique intersection.")

        x_x = (C1 * B2 - B1 * C2) / det
        x_y = (A1 * C2 - C1 * A2) / det

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
