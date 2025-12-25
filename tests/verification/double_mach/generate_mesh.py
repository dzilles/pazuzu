import gmsh
import sys
import os

def generate_mesh(filename):
    gmsh.initialize()
    gmsh.model.add("double_mach")

    # Domain [0, 4] x [0, 1]
    # Split bottom boundary at x = 1/6
    x_split = 1.0 / 6.0
    
    # Points
    p1 = gmsh.model.geo.addPoint(0, 0, 0)
    p2 = gmsh.model.geo.addPoint(x_split, 0, 0)
    p3 = gmsh.model.geo.addPoint(4, 0, 0)
    p4 = gmsh.model.geo.addPoint(4, 1, 0)
    p5 = gmsh.model.geo.addPoint(0, 1, 0)

    # Lines
    l1 = gmsh.model.geo.addLine(p1, p2) # Bottom left (Exact)
    l2 = gmsh.model.geo.addLine(p2, p3) # Bottom right (Wall)
    l3 = gmsh.model.geo.addLine(p3, p4) # Right (Exact)
    l4 = gmsh.model.geo.addLine(p4, p5) # Top (Exact)
    l5 = gmsh.model.geo.addLine(p5, p1) # Left (Exact)

    # Curve Loop
    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5])
    
    # Plane Surface
    s = gmsh.model.geo.addPlaneSurface([cl])

    # Transfinite mesh
    # Nx = 400, Ny = 100 for a decent resolution (4:1 aspect ratio)
    nx_total = 400
    ny = 100
    
    nx_left = int(nx_total * (x_split / 4.0))
    if nx_left < 2: nx_left = 2
    nx_right = nx_total - nx_left + 1
    
    gmsh.model.geo.mesh.setTransfiniteCurve(l1, nx_left)
    gmsh.model.geo.mesh.setTransfiniteCurve(l2, nx_right)
    gmsh.model.geo.mesh.setTransfiniteCurve(l3, ny)
    gmsh.model.geo.mesh.setTransfiniteCurve(l4, nx_total)
    gmsh.model.geo.mesh.setTransfiniteCurve(l5, ny)
    
    gmsh.model.geo.mesh.setTransfiniteSurface(s, "Left", [p1, p3, p4, p5])
    gmsh.model.geo.mesh.setRecombine(2, s)

    gmsh.model.geo.synchronize()
    
    # Physical Groups
    # Use IDs from src/kernels/boundary_conditions.py
    # BC_WALL = 1
    # BC_DOUBLE_MACH_EXACT = 7
    
    gmsh.model.addPhysicalGroup(1, [l1, l3, l4, l5], 7, "Exact")
    gmsh.model.addPhysicalGroup(1, [l2], 1, "Wall")
    gmsh.model.addPhysicalGroup(2, [s], 10, "Domain")
    
    gmsh.model.mesh.generate(2)
    
    output_dir = os.path.dirname(filename)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    gmsh.write(filename)
    gmsh.finalize()

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, "output", "double_mach.msh")
    generate_mesh(output_path)
