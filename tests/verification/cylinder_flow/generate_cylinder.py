import gmsh
import sys
import os

def generate_mesh(filename):
    gmsh.initialize()
    gmsh.model.add("cylinder")

    lc = 0.5
    lc_cyl = 0.1 # Finer near cylinder

    # Domain Points
    p1 = gmsh.model.geo.addPoint(-5, -5, 0, lc)
    p2 = gmsh.model.geo.addPoint(10, -5, 0, lc)
    p3 = gmsh.model.geo.addPoint(10, 5, 0, lc)
    p4 = gmsh.model.geo.addPoint(-5, 5, 0, lc)

    # Domain Lines
    l1 = gmsh.model.geo.addLine(p1, p2) # Bottom
    l2 = gmsh.model.geo.addLine(p2, p3) # Outlet
    l3 = gmsh.model.geo.addLine(p3, p4) # Top
    l4 = gmsh.model.geo.addLine(p4, p1) # Inlet

    # Cylinder (Radius 0.5)
    R = 0.5
    c1 = gmsh.model.geo.addPoint(R, 0, 0, lc_cyl)
    c2 = gmsh.model.geo.addPoint(0, R, 0, lc_cyl)
    c3 = gmsh.model.geo.addPoint(-R, 0, 0, lc_cyl)
    c4 = gmsh.model.geo.addPoint(0, -R, 0, lc_cyl)
    center = gmsh.model.geo.addPoint(0, 0, 0, lc_cyl)
    
    arc1 = gmsh.model.geo.addCircleArc(c1, center, c2)
    arc2 = gmsh.model.geo.addCircleArc(c2, center, c3)
    arc3 = gmsh.model.geo.addCircleArc(c3, center, c4)
    arc4 = gmsh.model.geo.addCircleArc(c4, center, c1)
    
    loop_cyl = gmsh.model.geo.addCurveLoop([arc1, arc2, arc3, arc4])
    loop_dom = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    
    s = gmsh.model.geo.addPlaneSurface([loop_dom, loop_cyl])

    gmsh.model.geo.synchronize()
    
    # Physical Groups
    gmsh.model.addPhysicalGroup(1, [l1], 1, "Bottom")
    gmsh.model.addPhysicalGroup(1, [l2], 2, "Outlet")
    gmsh.model.addPhysicalGroup(1, [l3], 3, "Top")
    gmsh.model.addPhysicalGroup(1, [l4], 4, "Inlet")
    gmsh.model.addPhysicalGroup(1, [arc1, arc2, arc3, arc4], 5, "Cylinder")
    gmsh.model.addPhysicalGroup(2, [s], 10, "Domain")
    
    # Mesh settings
    # To get good quads near cylinder, we might need a Boundary Layer or just basic recombination
    gmsh.model.mesh.setRecombine(2, s) # Quads
    gmsh.model.mesh.generate(2)
    
    # Ensure output dir exists
    output_dir = os.path.dirname(filename)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    gmsh.write(filename)
    gmsh.finalize()

if __name__ == "__main__":
    # Default to saving in 'output' subdirectory relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, "output", "cylinder.msh")
    generate_mesh(output_path)
