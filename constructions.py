## Build the problem in Python

from matplotlib import pyplot as plt
from typing import Dict, Tuple, Set
import numpy as np
from relations import *


def plot_point(ax, points: Dict[str, Tuple[float, float]]):
    for point in points.keys():
        ax.plot(points[point][0], points[point][1], 'o', label=point)

def plot_lines_from_eq(ax, lines: Dict[str, Tuple[float, float, float]]):
    for line in lines.keys():
        x_lim = ax.get_xlim()
        x_vals = np.array(x_lim)
        y_vals = (-lines[line][0] * x_vals - lines[line][2]) / lines[line][1]
        ax.plot(x_vals, y_vals)

def plot_circles_from_eq(ax, circles: Dict[str, Tuple[float, float, float]]):
    for circle in circles.keys():
        theta = np.linspace(0, 2 * np.pi, 200) 
        x_coords = circles[circle][0] + circles[circle][2] * np.cos(theta)
        y_coords = circles[circle][1] + circles[circle][2] * np.sin(theta)
        ax.plot(x_coords, y_coords)

def plot_midpoint(ax, name: str, a: Point, b: Point):
    new_x = (b.x - a.x)/2
    new_y = (b.y - b.x)/2
    ax.plot(new_x, new_y, 'o', label='name')

def plot_newline(ax, a: Point, b: Point):
    ax.plot([a.x, b.x], [a.y, b.y])

def plot_anglebisector(ax, a: Point, b: Point, c: Point):
    A = np.array([a.x, a.y])
    B = np.array([b.x, b.y])
    C = np.array([c.x, c.y])
    BA = A - B
    BC = C - B
    BA_norm = BA / np.linalg.norm(BA)
    BC_norm = BC / np.linalg.norm(BC)
    bisector_dir = BA_norm + BC_norm
    bisector_dir = bisector_dir / np.linalg.norm(bisector_dir)
    start_point = B - bisector_dir * 100
    end_point = B + bisector_dir * 100
    plt.plot([start_point[0], end_point[0]], [start_point[1], end_point[1]])

def plot_perp(ax, a: Point, b: Point, pt_from: Point):
    A = np.array([a.x, a.y])
    B = np.array([b.x, b.y])
    AB = A - B
    perp = np.array([-AB[1], AB[0]])
    perp = perp / np.linalg.norm(perp)
    start_point = pt_from - perp * 100
    end_point = pt_from + perp * 100
    plt.plot([start_point[0], end_point[0]], [start_point[1], end_point[1]])