# The MIT License (MIT)
#
# Copyright (c) 2022 Daniel Zilles
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

import gmsh
import numpy

## Class containing all function and objects necessary for building the mesh:
#
#  Example
#  ~~~~~~~~~~~~~~~~~~~~~
#  ~~~~~~~~~~~~~~~~~~~~~
#
class Mesh:

    ## Initialization of the class.
    #  @param self The object pointer.
    #
    def __init__(self, path_to_file):

        ## Coordinates of the mesh
        #self.coordinates = numpy.array([], dtype=float);

        ## Indices of nodes yielding cells
        #self.cell_node_ids = numpy.array([], dtype=int);

        ## Indices of nodes yielding faces
        #self.face_node_ids = numpy.array([], dtype=int);

        ## Indices of neighbor nodes
        #self.neighbor_ids = numpy.array([], dtype=int)

        ## Indices of adjacent faces
        #self.adjacent_face_ids = numpy.array([], dtype=int)

        ## Load the data from a gmsh file
        self.load_gmsh(path_to_file)

    ## Set new data
    #  @param self The object pointer.
    #  @param data Numpy matrix with the data values.
    #
    #  This function loads the mesh data of a gmsh file.
    #
    def load_gmsh(self, path_to_file):

        # Open gmsh file
        gmsh.initialize();
        gmsh.open(path_to_file);
        tmp1 = gmsh.model.mesh.getNodes();
        tmp2 = gmsh.model.mesh.getNodesByElementType(2);

        #print(tmp1[1].reshape(-1,3))
         
        gmsh.finalize();

        self.number_nodes = tmp1[1].reshape(-1,3).shape[0];
        self.number_cells = tmp2[0].reshape(-1,3).shape[0];
        self.number_faces_per_cell =  3;

        # Store the coordinates in internal structure
        self.coordinates = numpy.zeros((self.number_nodes, 3), dtype=int);
 
        #self.coordinates[:,0] = tmp1[1].reshape(-1,3)[:,0];
        #self.coordinates[:,1] = tmp1[1].reshape(-1,3)[:,1];
        #self.coordinates[:,2] = tmp1[1].reshape(-1,3)[:,2];
 
        self.coordinates = tmp1[1].reshape(-1,3)

        #print(self.coordinates)

        # Store the nodes in internal structure
        self.cell_node_ids = numpy.zeros((self.number_cells, 3), dtype=int);

        #self.cell_node_ids[:,0] = tmp2[0].reshape(-1,3)[:,0];
        #self.cell_node_ids[:,1] = tmp2[0].reshape(-1,3)[:,1];
        #self.cell_node_ids[:,2] = tmp2[0].reshape(-1,3)[:,2];

        self.cell_node_ids = tmp2[0].reshape(-1,3)

        # Build faces using the cell node ids
        self.face_node_ids = numpy.zeros((3*self.number_cells, 2), dtype=int);

        self.face_node_ids[0:self.number_cells,:] = self.cell_node_ids[:,0:2];
        self.face_node_ids[self.number_cells:2*self.number_cells,:] = self.cell_node_ids[:,1:3];
        self.face_node_ids[2*self.number_cells:3*self.number_cells,:] = numpy.vstack((self.cell_node_ids[:,2], self.cell_node_ids[:,0])).transpose();

        self.face_node_ids.sort(axis=1);

        element_to_element = numpy.linspace(1, self.number_cells, self.number_cells).reshape(1,-1)*numpy.ones(self.number_faces_per_cell).reshape(-1,1)

        element_to_face = numpy.ones(self.number_cells).reshape(1,-1)*numpy.linspace(1, self.number_faces_per_cell, self.number_faces_per_cell).reshape(-1,1)

        ids = self.face_node_ids[:,0]*self.number_nodes + self.face_node_ids[:,1] +1;

        spNodeToNode = numpy.zeros((self.number_cells*self.number_faces_per_cell, 4));

        spNodeToNode[:,0] = ids;
        spNodeToNode[:,1] = numpy.linspace(1, self.number_faces_per_cell*self.number_cells, self.number_faces_per_cell*self.number_cells);
        spNodeToNode[:,2] = element_to_element.flatten();
        spNodeToNode[:,3] = element_to_face.flatten();


        spNodeToNode = spNodeToNode[spNodeToNode[:, 3].argsort()]
        spNodeToNode = spNodeToNode[spNodeToNode[:, 2].argsort(kind='mergesort')]
        spNodeToNode = spNodeToNode[spNodeToNode[:, 1].argsort(kind='mergesort')]
        spNodeToNode = spNodeToNode[spNodeToNode[:, 0].argsort(kind='mergesort')]

        tmp = spNodeToNode[:,0] - numpy.roll(spNodeToNode[:,0], -1) 
   
        indices = numpy.where(tmp == 0)[0]

        matchL = numpy.vstack((spNodeToNode[indices], spNodeToNode[indices+1])).astype(int)
        matchR = numpy.vstack((spNodeToNode[indices+1], spNodeToNode[indices])).astype(int)

        element_to_element = element_to_element.transpose()

        for i in range(0, 2*indices.size):
 
          element_to_element[matchL[i, 2]-1, matchL[i, 3]-1] = matchR[i,2];
          element_to_face.reshape(-1,1)[matchL[i, 1]-1] = matchR[i,3];

        self.neighbor_ids = (element_to_element - 1).astype(int);
        self.adjacent_face_ids = (element_to_face - 1).astype(int);
        self.cell_node_ids = (self.cell_node_ids - 1).astype(int);
        self.face_node_ids = (self.face_node_ids - 1).astype(int);
