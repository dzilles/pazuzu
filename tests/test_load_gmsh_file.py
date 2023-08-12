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

