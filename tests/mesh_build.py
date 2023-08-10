# Import modules:
import gmsh
import sys

# Initialize gmsh:
gmsh.initialize()

# cube points:
lc = 1e-1
point1 = gmsh.model.geo.add_point(0, 0, 0, lc)
point2 = gmsh.model.geo.add_point(0.3, 0, 0, lc)
point3 = gmsh.model.geo.add_point(0.3, 0.3, 0, lc)
point4 = gmsh.model.geo.add_point(0, 0.3, 0, lc)

# Edge of cube:
line1 = gmsh.model.geo.add_line(point1, point2)
line2 = gmsh.model.geo.add_line(point2, point3)
line3 = gmsh.model.geo.add_line(point3, point4)
line4 = gmsh.model.geo.add_line(point4, point1)

# faces of cube:
face1 = gmsh.model.geo.add_curve_loop([line1, line2, line3, line4])

# surfaces of cube:
gmsh.model.geo.add_plane_surface([face1])

# Create the relevant Gmsh data structures
# from Gmsh model.
gmsh.model.geo.synchronize()

# Generate mesh:
gmsh.model.mesh.generate()

# Write mesh data:
gmsh.write("GFG.msh")

# Creates graphical user interface
if 'close' not in sys.argv:
	gmsh.fltk.run()

# It finalize the Gmsh API
gmsh.finalize()

