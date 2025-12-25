import gmsh
import sys
import os

def generate_mesh(filename):
    gmsh.initialize()
    gmsh.model.add("sod_shock_tube")

    # Domain [0, 1.0] x [0, 0.1]
    # Resolution: element size 0.01
    
    # Points
    p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, 0.01)
    p2 = gmsh.model.geo.addPoint(1.0, 0.0, 0.0, 0.01)
    p3 = gmsh.model.geo.addPoint(1.0, 0.1, 0.0, 0.01)
    p4 = gmsh.model.geo.addPoint(0.0, 0.1, 0.0, 0.01)

    # Lines
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    # Curve Loop
    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    
    # Plane Surface
    s = gmsh.model.geo.addPlaneSurface([cl])

    # Transfinite mesh to get perfect quads
    # 100 elements in X -> 101 nodes
    # 10 elements in Y -> 11 nodes
    gmsh.model.geo.mesh.setTransfiniteCurve(l1, 101)
    gmsh.model.geo.mesh.setTransfiniteCurve(l2, 11)
    gmsh.model.geo.mesh.setTransfiniteCurve(l3, 101)
    gmsh.model.geo.mesh.setTransfiniteCurve(l4, 11)
    
    gmsh.model.geo.mesh.setTransfiniteSurface(s)
    gmsh.model.geo.mesh.setRecombine(2, s) # Recombine triangles into quads

    gmsh.model.geo.synchronize()
    
    # Add Physical Groups
    gmsh.model.addPhysicalGroup(1, [l1], 1, "Bottom")
    gmsh.model.addPhysicalGroup(1, [l2], 2, "Right")
    gmsh.model.addPhysicalGroup(1, [l3], 3, "Top")
    gmsh.model.addPhysicalGroup(1, [l4], 4, "Left")
    gmsh.model.addPhysicalGroup(2, [s], 10, "Domain")
    
    gmsh.model.mesh.generate(2)
    
    # Ensure output dir exists
    output_dir = os.path.dirname(filename)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    gmsh.write(filename)
    gmsh.finalize()

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, "output", "sod.msh")
    generate_mesh(output_path)
