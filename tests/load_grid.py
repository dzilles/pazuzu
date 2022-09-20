import gmsh
import numpy
from pazuzu import *

gmsh.initialize();

gmsh.open("vertex_2D.msh");

a = gmsh.model.mesh.getNodesByElementType(2);

print(a[0]);

mesh = Mesh(30, 30);

print(mesh.node_count);

b = a[1].reshape(-1, 3);

mesh.set_coordinates(b, a[0].reshape(-1, 3));


gmsh.finalize();

