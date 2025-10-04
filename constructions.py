from matplotlib import pyplot as plt
from typing import Dict, Tuple
import numpy as np
from relations import *

plt.style.use("seaborn-v0_8-whitegrid")  # clean white background


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
        ax.text(x, y, name, fontsize=9, color='black', ha='center', va='center')


def plot_midpoint(ax, name: str, a: Point, b: Point):
    new_x = (a.x + b.x) / 2
    new_y = (a.y + b.y) / 2
    ax.plot(new_x, new_y, 'o', color='black', markersize=6, zorder=5)
    ax.text(new_x + 0.05, new_y + 0.05, name, fontsize=9, color='black', weight='bold')


def plot_newline(ax, a: Point, b: Point):
    ax.plot([a.x, b.x], [a.y, b.y], color='black', lw=1.0, alpha=0.9, linestyle='-')


def plot_anglebisector(ax, a: Point, b: Point, c: Point):
    A, B, C = np.array([a.x, a.y]), np.array([b.x, b.y]), np.array([c.x, c.y])
    BA, BC = A - B, C - B
    BA_norm = BA / np.linalg.norm(BA)
    BC_norm = BC / np.linalg.norm(BC)
    bisector_dir = (BA_norm + BC_norm) / np.linalg.norm(BA_norm + BC_norm)
    start_point = B - bisector_dir * 100
    end_point = B + bisector_dir * 100
    ax.plot([start_point[0], end_point[0]], [start_point[1], end_point[1]],
            color='black', lw=1.0, linestyle='--', alpha=0.9)


def plot_perp(ax, a: Point, b: Point, pt_from: Point):
    A, B = np.array([a.x, a.y]), np.array([b.x, b.y])
    AB = A - B
    perp = np.array([-AB[1], AB[0]])
    perp /= np.linalg.norm(perp)
    start_point = np.array([pt_from.x, pt_from.y]) - perp * 100
    end_point = np.array([pt_from.x, pt_from.y]) + perp * 100
    ax.plot([start_point[0], end_point[0]], [start_point[1], end_point[1]],
            color='black', lw=1.0, linestyle=':', alpha=0.9)


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
