import gmsh
import sys
import os

def generate_mesh(filename):
    gmsh.initialize()
    gmsh.model.add("pulse")

    # Domain [-1, 1] x [-1, 1]
    L = 1.0
    
    # Points
    p1 = gmsh.model.geo.addPoint(-L, -L, 0)
    p2 = gmsh.model.geo.addPoint(L, -L, 0)
    p3 = gmsh.model.geo.addPoint(L, L, 0)
    p4 = gmsh.model.geo.addPoint(-L, L, 0)

    # Lines
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    # Curve Loop
    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    
    # Plane Surface
    s = gmsh.model.geo.addPlaneSurface([cl])

    # Transfinite mesh
    # dx ~ 0.05 -> 40 elements per side (2.0 / 40 = 0.05)
    n = 41
    gmsh.model.geo.mesh.setTransfiniteCurve(l1, n)
    gmsh.model.geo.mesh.setTransfiniteCurve(l2, n)
    gmsh.model.geo.mesh.setTransfiniteCurve(l3, n)
    gmsh.model.geo.mesh.setTransfiniteCurve(l4, n)
    
    gmsh.model.geo.mesh.setTransfiniteSurface(s)
    gmsh.model.geo.mesh.setRecombine(2, s)

    gmsh.model.geo.synchronize()
    
    # Physical Groups
    gmsh.model.addPhysicalGroup(1, [l1], 1, "Bottom")
    gmsh.model.addPhysicalGroup(1, [l2], 2, "Right")
    gmsh.model.addPhysicalGroup(1, [l3], 3, "Top")
    gmsh.model.addPhysicalGroup(1, [l4], 4, "Left")
    gmsh.model.addPhysicalGroup(2, [s], 10, "Domain")
    
    gmsh.model.mesh.generate(2)
    
    output_dir = os.path.dirname(filename)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    gmsh.write(filename)
    gmsh.finalize()

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, "output", "pulse.msh")
    generate_mesh(output_path)
