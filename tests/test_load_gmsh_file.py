from pazuzu import mesh
import numpy
import pytest

#from pazuzu import *
#from pazuzu import visualization as vis

build_test_cases = True

def test_load_2D_gmsh_coordinates():

    msh = mesh.Mesh("vertex_2D.msh");

    if (build_test_cases):
        numpy.save("test_data/test_load_gmsh_coordinates", msh.coordinates)

    test_data = numpy.load("test_data/test_load_gmsh_coordinates.npy")

    assert (numpy.array_equal(msh.coordinates, test_data))


          



#    assert True



#import numpy
#from pazuzu import *
#from pazuzu import visualization as vis
#import sys
#numpy.set_printoptions(threshold=sys.maxsize)




#v = vis.Visualization();
#v.initialize_scalar_bar();

#v.build_grid(msh.coordinates[:,0], msh.coordinates[:,1], msh.cell_node_ids)

#data = numpy.ones((1868,1)).astype(float)*0.5

#aa = numpy.vstack((numpy.arange(1868), msh.cell_node_ids.transpose()))

#print(msh.cell_node_ids)
#print(msh.coordinates)

#for i in range(0,42):
#    data = numpy.ones((51868)).astype(float)*1.0
#    print(msh.cell_node_ids[i])
#    print(msh.coordinates[10])
#    print(msh.coordinates[11])
#    print(msh.coordinates[13])
#    print(i, msh.neighbor_ids[i,0], msh.neighbor_ids[i,1], msh.neighbor_ids[i,2])
#    data[msh.neighbor_ids[i,0]] = 0.2
#    data[msh.neighbor_ids[i,1]] = 0.2
#    data[msh.neighbor_ids[i,2]] = 0.2
#    data[i] = 0.5
#    v.plot_data(data, [0,1])
#    input("Press Enter to continue...")

#data[2] = 0.2;
#data[480] = 0.;
#data[481] = 0.0;
#data[938] = 0.0;


#v.plot_grid(msh.coordinates[:,0], msh.coordinates[:,1], msh.cell_node_ids)

#while(True):
#    v.Render()
