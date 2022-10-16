import gmsh
import numpy
from pazuzu import *
from pazuzu import visualization as vis
import sys
numpy.set_printoptions(threshold=sys.maxsize)


gmsh.initialize();

gmsh.open("vertex_2D.msh");

a = gmsh.model.mesh.getNodesByElementType(2);

#print(a[0]);

#mesh = Mesh(30, 30);

#print(mesh.node_count);

#b = a[1].reshape(-1, 3);

#mesh.set_coordinates(b, a[0].reshape(-1, 3));

x = a[1].reshape(-1, 3)[:,0];
y = a[1].reshape(-1, 3)[:,1];

v = vis.Visualization();
#v.initialize_scalar_bar();

#v.initialize_data_plot(x, y)

#v.plot_grid(x, y);

#while(True):
#v.Render()

v = numpy.arange(0, a[0].reshape(-1,3).shape[0], 1).reshape(-1,1)
vv = a[0].reshape(-1,3)

#vv.sort(axis=1)

#print(vv)

number_cells = a[0].reshape(-1,3).shape[0];
number_faces = 3;
number_nodes = x.shape[0]

#boundary_condition_ids = numpy.zeros(number_cells, number_faces);

#for f in range(0, number_faces):

#  for c in range(0, number_cells):

#    if()


#print(vv)
tmp1 = numpy.column_stack((vv[:,0], vv[:,1]))
tmp2 = numpy.column_stack((vv[:,1], vv[:,2]))
tmp3 = numpy.column_stack((vv[:,2], vv[:,0]))
 
fnodes = numpy.row_stack((tmp1, tmp2, tmp3)) 

fnodes.sort(axis=1)
   

EToE = numpy.arange(0, number_cells, 1).reshape(-1,1)*numpy.ones(number_faces).reshape(1, number_faces).astype(int)

EToF = numpy.ones(number_cells).reshape(-1,1)*numpy.arange(0, number_faces).reshape(1, number_faces)

ident = fnodes[:,0]*number_nodes + fnodes[:,1]

ident2 = numpy.arange(0,number_faces*number_cells)

spNodeToNode = (numpy.column_stack((ident.reshape(-1,1), ident2.reshape(-1,1), EToE.reshape(-1,1), EToF.reshape(-1,1)))).astype(int)

spNodeToNode = spNodeToNode[numpy.lexsort(spNodeToNode.T[::-1])]

indices = (numpy.where(numpy.equal(spNodeToNode[0:-2, 0], spNodeToNode[1:-1,0]))).reshape(-1,1)

inices2 = indices + numpy.ones(indices.shape())


matchL = spNodeToNode[indices,:] 
#matchR = spNodeToNode[indices+1,:]

print(matchL)

gmsh.finalize();

