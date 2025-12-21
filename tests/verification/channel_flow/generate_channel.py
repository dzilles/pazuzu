import gmsh
import sys
import os

def generate_mesh(filename):
    gmsh.initialize()
    gmsh.model.add("channel")

    L = 10.0
    H = 2.0
    lc = 0.5
    
    # Points
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(L, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(L, H, 0, lc)
    p4 = gmsh.model.geo.addPoint(0, H, 0, lc)

    # Lines
    l1 = gmsh.model.geo.addLine(p1, p2) # Bottom
    l2 = gmsh.model.geo.addLine(p2, p3) # Outlet
    l3 = gmsh.model.geo.addLine(p3, p4) # Top
    l4 = gmsh.model.geo.addLine(p4, p1) # Inlet

    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    s = gmsh.model.geo.addPlaneSurface([cl])

    # Structured grid
    gmsh.model.geo.mesh.setTransfiniteCurve(l1, 21)
    gmsh.model.geo.mesh.setTransfiniteCurve(l3, 21)
    gmsh.model.geo.mesh.setTransfiniteCurve(l2, 6)
    gmsh.model.geo.mesh.setTransfiniteCurve(l4, 6)
    
    gmsh.model.geo.mesh.setTransfiniteSurface(s)
    gmsh.model.geo.mesh.setRecombine(2, s)

    gmsh.model.geo.synchronize()
    
    # Physical Groups
    gmsh.model.addPhysicalGroup(1, [l1], 1, "WallBottom")
    gmsh.model.addPhysicalGroup(1, [l2], 2, "Outlet")
    gmsh.model.addPhysicalGroup(1, [l3], 3, "WallTop")
    gmsh.model.addPhysicalGroup(1, [l4], 4, "Inlet")
    gmsh.model.addPhysicalGroup(2, [s], 10, "Domain")
    
    gmsh.model.mesh.generate(2)
    gmsh.write(filename)
    gmsh.finalize()

if __name__ == "__main__":
    generate_mesh("channel.msh")
