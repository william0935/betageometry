import os
import zipfile
from typing import Dict, Tuple, Set
import xml.etree.ElementTree as ET
import math


def invert_2x2(m, tol=1e-12):
    (a, b), (c, d) = m
    det = a * d - b * c
    if abs(det) < tol:
        return None
    return [[d / det, -b / det],
            [-c / det, a / det]]


def parse_conic_as_circle(A, B, C, D, E, F, tol=1e-9):
    """
    Interpret the conic A*x^2 + B*x*y + C*y^2 + D*x + E*y + F=0
    as a circle (possibly rotated), returning (cx, cy, r) or None if it fails.

    Conditions for unrotated circle (approx):
     - A ~ C (within tol)
     - B^2 < 4*A*C
     - radius^2 > 0
    """
    if abs(A - C) > tol:
        return None
    if abs(A) < tol:
        return None
    if B * B >= 4 * A * C:
        return None

    M = [[A, B / 2],
         [B / 2, C]]
    v = [D, E]
    invM = invert_2x2(M, tol=tol)
    if not invM:
        return None

    # center = -0.5 * invM * v
    cx = -0.5 * (invM[0][0] * v[0] + invM[0][1] * v[1])
    cy = -0.5 * (invM[1][0] * v[0] + invM[1][1] * v[1])

    # radius^2 = (cx,cy)*M*(cx,cy)^T - F
    Mx = A * cx + (B / 2) * cy
    My = (B / 2) * cx + C * cy
    r_sq = (cx * Mx + cy * My) - F
    if r_sq <= 0:
        return None

    r = math.sqrt(r_sq)
    return (cx, cy, r)


def circle_from_3_points(p1, p2, p3, tol=1e-12):
    """
    Given 3 distinct points (x1,y1), (x2,y2), (x3,y3),
    return the circle center (cx,cy) and radius r.
    If collinear or degenerate, return None.
    """
    (x1, y1) = p1
    (x2, y2) = p2
    (x3, y3) = p3

    d = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    if abs(d) < tol:
        return None  # Points nearly collinear

    x1_sq = x1 * x1 + y1 * y1
    x2_sq = x2 * x2 + y2 * y2
    x3_sq = x3 * x3 + y3 * y3

    cx = (x1_sq * (y2 - y3) + x2_sq * (y3 - y1) + x3_sq * (y1 - y2)) / d
    cy = (x1_sq * (x3 - x2) + x2_sq * (x1 - x3) + x3_sq * (x2 - x1)) / d
    r = math.dist((cx, cy), (x1, y1))
    return (cx, cy, r)

def parse_ggb_file(ggb_filename):
    """
    Parses a GeoGebra (.ggb) file to extract geometric data.
    Returns dictionaries for points, lines, segments, and circles.
    """

    def circle_from_center_and_point(center, point):
        """Helper to calculate circle data from center and point."""
        cx, cy = center
        px, py = point
        radius = math.dist((cx, cy), (px, py))
        return cx, cy, radius

    # Prepare extraction folder
    project_dir = os.path.dirname(os.path.abspath(__file__))
    extract_folder = os.path.join(project_dir, "temp_ggb_extraction")
    os.makedirs(extract_folder, exist_ok=True)

    # Unzip the .ggb file
    try:
        with zipfile.ZipFile(ggb_filename, 'r') as zf:
            zf.extractall(extract_folder)
    except zipfile.BadZipFile:
        raise ValueError("Invalid .ggb file (not a valid ZIP).")

    # Parse geogebra.xml
    xml_path = os.path.join(extract_folder, "geogebra.xml")
    if not os.path.exists(xml_path):
        raise FileNotFoundError("No geogebra.xml found inside the .ggb archive.")

    print(xml_path)

    tree = ET.parse(xml_path)
    root = tree.getroot()
    construction = root.find('construction')
    if construction is None:
        return {}, {}, {}, {}

    # A) Collect Points
    points_dict = {}
    for elem in construction.findall('element'):
        if elem.get('type') == 'point':
            label = elem.get('label')
            coords = elem.find('coords')
            if label and coords is not None:
                try:
                    x_val = float(coords.get("x", "0"))
                    y_val = float(coords.get("y", "0"))
                    z_val = float(coords.get("z", "0"))
                    points_dict[label] = (x_val/z_val, y_val/z_val)
                except ValueError:
                    pass

    # B) Lines & Segments
    lines_dict = {}
    segments_dict = {}
    for elem in construction.findall('element'):
        etype = elem.get('type')
        label = elem.get('label')
        if etype in ("line", "segment") and label:
            cnode = elem.find('coords')
            if cnode is not None:
                try:
                    A_ = float(cnode.get("x", "0"))
                    B_ = float(cnode.get("y", "0"))
                    C_ = float(cnode.get("z", "0"))
                    if etype == "line":
                        lines_dict[label] = (A_, B_, C_)
                    else:
                        segments_dict[label] = (A_, B_, C_)
                except ValueError:
                    pass

    # C) Conics → Circles
    conic_circles = {}
    for elem in construction.findall('element'):
        if elem.get('type') == 'conic':
            label = elem.get('label')
            mat = elem.find('matrix')
            if label and mat is not None:
                try:
                    A0 = float(mat.get("A0", "0"))
                    A1 = float(mat.get("A1", "0"))
                    A2 = float(mat.get("A2", "0"))
                    A3 = float(mat.get("A3", "0"))
                    A4 = float(mat.get("A4", "0"))
                    A5 = float(mat.get("A5", "0"))
                except ValueError:
                    continue
                circle_data = parse_conic_as_circle(A0, A1, A2, A3, A4, A5)
                if circle_data:
                    conic_circles[label] = circle_data

    # D) Circles from <command name="Circle">
    command_circles = {}
    for cmd in construction.findall('command'):
        if cmd.get('name') == 'Circle':
            inputs = cmd.find('input')
            outputs = cmd.find('output')
            if inputs is None or outputs is None:
                continue
            in_labels = [v for k, v in sorted(inputs.items()) if k.startswith('a')]
            out_labels = [v for k, v in sorted(outputs.items()) if k.startswith('a')]

            if len(out_labels) == 1:
                circ_label = out_labels[0]
            else:
                continue

            if len(in_labels) == 3:  # Circle defined by 3 points
                pA, pB, pC = in_labels
                if all(lbl in points_dict for lbl in (pA, pB, pC)):
                    coordsA = points_dict[pA]
                    coordsB = points_dict[pB]
                    coordsC = points_dict[pC]
                    circle_data = circle_from_3_points(coordsA, coordsB, coordsC)
                    if circle_data:
                        command_circles[circ_label] = circle_data

            elif len(in_labels) == 2:  # Circle defined by center and point
                center_label, point_label = in_labels
                if center_label in points_dict and point_label in points_dict:
                    center_coords = points_dict[center_label]
                    point_coords = points_dict[point_label]
                    circle_data = circle_from_center_and_point(center_coords, point_coords)
                    if circle_data:
                        command_circles[circ_label] = circle_data

    # Merge circles from conics and commands
    final_circles = {**conic_circles, **command_circles}

    return points_dict, lines_dict, segments_dict, final_circles

def rename_circle_dictionary(
    circles: Dict[str, Tuple[float, float, float]]
) -> Dict[str, Tuple[float, float, float]]:
    """
    Renames the keys in the circle dictionary to sequential labels 'Circle1', 'Circle2', etc.

    Args:
        circles (Dict[str, Tuple[float, float, float]]): A dictionary with circle data where
            keys are labels and values are tuples of the form (cx, cy, r).

    Returns:
        Dict[str, Tuple[float, float, float]]: A new dictionary with renamed keys.

    Raises:
        TypeError: If the input is not a dictionary.
    """
    if not isinstance(circles, dict):
        raise TypeError("Input must be a dictionary.")

    renamed_circles = {}
    for idx, (key, value) in enumerate(circles.items(), start=1):
        renamed_circles[f"Circle{idx}"] = value

    return renamed_circles

def combine_lines_and_segments_dictionaries(
    dict1: Dict[str, Tuple[float, float, float]],
    dict2: Dict[str, Tuple[float, float, float]]
) -> Dict[str, Tuple[float, float, float]]:
    """
    Combines two dictionaries and renames the keys sequentially as 'L1', 'L2', 'L3', etc.

    Args:
        dict1 (Dict[str, Tuple[float, float, float]]): The first dictionary.
        dict2 (Dict[str, Tuple[float, float, float]]): The second dictionary.

    Returns:
        Dict[str, Tuple[float, float, float]]: A new dictionary with keys renamed sequentially.

    Raises:
        TypeError: If either dict1 or dict2 is not a dictionary.
    """
    if not isinstance(dict1, dict) or not isinstance(dict2, dict):
        raise TypeError("Both inputs must be dictionaries.")

    combined_dict = {}
    all_items = list(dict1.items()) + list(dict2.items())

    for idx, (_, value) in enumerate(all_items, start=1):
        combined_dict[f"L{idx}"] = value

    return combined_dict


def rename_circle_dictionary(
    circles: Dict[str, Tuple[float, float, float]]
) -> Dict[str, Tuple[float, float, float]]:
    """
    Renames the keys in the circle dictionary to sequential labels 'Circle1', 'Circle2', etc.

    Args:
        circles (Dict[str, Tuple[float, float, float]]): A dictionary with circle data where
            keys are labels and values are tuples of the form (cx, cy, r).

    Returns:
        Dict[str, Tuple[float, float, float]]: A new dictionary with renamed keys.

    Raises:
        TypeError: If the input is not a dictionary.
    """
    if not isinstance(circles, dict):
        raise TypeError("Input must be a dictionary.")

    renamed_circles = {}
    for idx, (key, value) in enumerate(circles.items(), start=1):
        renamed_circles[f"Circle{idx}"] = value

    return renamed_circles

def combine_lines_and_segments_dictionaries(
    dict1: Dict[str, Tuple[float, float, float]],
    dict2: Dict[str, Tuple[float, float, float]]
) -> Dict[str, Tuple[float, float, float]]:
    """
    Combines two dictionaries and renames the keys sequentially as 'L1', 'L2', 'L3', etc.

    Args:
        dict1 (Dict[str, Tuple[float, float, float]]): The first dictionary.
        dict2 (Dict[str, Tuple[float, float, float]]): The second dictionary.

    Returns:
        Dict[str, Tuple[float, float, float]]: A new dictionary with keys renamed sequentially.

    Raises:
        TypeError: If either dict1 or dict2 is not a dictionary.
    """
    if not isinstance(dict1, dict) or not isinstance(dict2, dict):
        raise TypeError("Both inputs must be dictionaries.")

    combined_dict = {}
    all_items = list(dict1.items()) + list(dict2.items())

    for idx, (_, value) in enumerate(all_items, start=1):
        combined_dict[f"L{idx}"] = value

    return combined_dict




def parse_picture(filename: str):
    # Build a path to the .ggb file inside "Geogebra Files"
    script_dir = os.path.dirname(os.path.abspath(__file__))
    ggb_path = os.path.join(script_dir, "Geogebra Files", filename)
    points, lines, segments, circles = parse_ggb_file(ggb_path)
    new_circles = rename_circle_dictionary(circles)
    new_lines = combine_lines_and_segments_dictionaries(lines,segments)
    return points, new_lines, new_circles

x,y,z= parse_picture('test.ggb')
print(x)
print(y)
print(z)