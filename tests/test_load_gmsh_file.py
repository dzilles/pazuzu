from pazuzu import mesh
import numpy
import pytest

#from pazuzu import *
#from pazuzu import visualization as vis

build_test_cases = False

def test_load_2D_gmsh_coordinates():

    msh = mesh.Mesh("vertex_2D.msh");

    if (build_test_cases):
        numpy.save("data/system_tests/test_load_gmsh_coordinates", msh.coordinates)
        numpy.save("data/system_tests/test_load_gmsh_cell_node_ids", msh.cell_node_ids)
        numpy.save("data/system_tests/test_load_gmsh_face_node_ids", msh.face_node_ids)
        numpy.save("data/system_tests/test_load_gmsh_neighbor_ids", msh.neighbor_ids)
        numpy.save("data/system_tests/test_load_gmsh_adjacent_face_ids", msh.adjacent_face_ids)

    test_data = numpy.load("data/system_tests/test_load_gmsh_coordinates.npy")
    assert (numpy.array_equal(msh.coordinates, test_data))

    test_data = numpy.load("data/system_tests/test_load_gmsh_cell_node_ids.npy")
    assert (numpy.array_equal(msh.cell_node_ids, test_data))

    test_data = numpy.load("data/system_tests/test_load_gmsh_face_node_ids.npy")
    assert (numpy.array_equal(msh.face_node_ids, test_data))

    test_data = numpy.load("data/system_tests/test_load_gmsh_neighbor_ids.npy")
    assert (numpy.array_equal(msh.neighbor_ids, test_data))

    test_data = numpy.load("data/system_tests/test_load_gmsh_adjacent_face_ids.npy")
    assert (numpy.array_equal(msh.adjacent_face_ids, test_data))
