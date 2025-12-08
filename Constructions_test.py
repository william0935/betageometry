from Constructions import Canva, Point
import numpy as np

def random_point_coords():
    return np.random.rand() * 10, np.random.rand() * 10

def create_sample_points():
    pts = []
    for name in ["A", "B", "C"]:
        x, y = random_point_coords()
        pts.append(Point(name, x, y))
    return pts

def test_add_point(canva):
    print("Testing add_point...")
    p = canva.add_point(1.0, 2.0)
    print(f"Added point: {p.name} at ({p.x}, {p.y})")

def test_midpoint(canva, pts):
    print("Testing midpoint...")
    p, rels = canva.midpoint(pts[0], pts[1])
    print(f"Midpoint: {p.name} at ({p.x}, {p.y}), relations: {rels}")

def test_mirror(canva, pts):
    print("Testing mirror...")
    p, rels = canva.mirror(pts[0], pts[1])
    print(f"Mirror: {p.name} at ({p.x}, {p.y}), relations: {rels}")

def test_anglebisector(canva, pts):
    print("Testing anglebisector...")
    p, rels = canva.anglebisector(pts[0], pts[1], pts[2])
    print(f"Angle bisector point: {p}, relations: {rels}")

def test_foot(canva, pts):
    print("Testing foot...")
    p, rels = canva.foot(pts[0], pts[1], pts[2])
    print(f"Foot point: {p}, relations: {rels}")

def test_circle(canva, pts):
    print("Testing circle...")
    p, rels = canva.circle(pts[0], pts[1], pts[2])
    print(f"Circle center: {p}, relations: {rels}")

def test_incenter(canva, pts):
    print("Testing incenter...")
    pts_list, rels = canva.incenter(pts[0], pts[1], pts[2])
    print(f"Incenter points: {pts_list}, relations: {rels}")

def test_incenter2(canva, pts):
    print("Testing incenter2...")
    pts_list, rels = canva.incenter2(pts[0], pts[1], pts[2])
    print(f"Incenter2 points: {pts_list}, relations: {rels}")

def test_plot(canva):
    print("Testing plot...")
    canva.plot()
    print("Plot displayed successfully.")

def run_all_tests():
    # Create sample points
    pts = create_sample_points()
    
    # Initialize Canva
    points_dict = {p.name: (p.x, p.y) for p in pts}
    canva = Canva(points=pts, points_dict=points_dict, lines={}, circles={})
    
    # Run all tests
    # test_add_point(canva)
    # test_midpoint(canva, pts)
    # test_mirror(canva, pts)
    # test_anglebisector(canva, pts)
    # test_foot(canva, pts)
    # test_circle(canva, pts)
    # test_incenter(canva, pts)
    # test_incenter2(canva, pts)
    test_plot(canva)

if __name__ == "__main__":
    run_all_tests()
