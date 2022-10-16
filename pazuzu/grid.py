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

import vtk
from vtk.util import numpy_support

## Class containing all function and objects necessary for building the grid:
#
#  Example
#  ~~~~~~~~~~~~~~~~~~~~~
#  # Import the class
#  from pazuzu.utility import visualization as vis
# 
#  # Construction of the visualization object
#  v = vis.Visualization()
#  # Initialize scalar bar
#  v.initialize_scalar_bar()
#  # Initialize the data plot, x and y are numpy matrices
#  # containing the coordinates corresponding to the data points
#  v.initialize_data_plot(x, y)
#
#  # Plot the data, data is numpy array with the same size as x and y
#  # The data is plotted for values between 0.8 and 1.2
#  v.update_data_plot(data, [0.8, 1.2])
#  # Write the render image to a png file
#  v.WriteImage("path/to/output/output.png")
#  ~~~~~~~~~~~~~~~~~~~~~
#
class Grid:

    ## Initialization of the class.
    #  @param self The object pointer.
    #
    def __init__(self):

        ## VTK object to specify the behavior of a rendering window.
        self.render_window = vtk.vtkRenderWindow()
        self.render_window.AddRenderer(self.renderer)
        self.render_window.SetSize(1920,1080)

    ## Set new data and render the image.
    #  @param self The object pointer.
    #  @param data Numpy matrix with the data values.
    #  @param scalar_range Range of the data. Should be a list with 2 entries or None.
    #
    #  This function plots the data.
    #  The scalar_range is contains the range of the couloring scheme.
    #  If the scalar_range is None, the minimum and maximum number in the data matrix are chosen.
    #
    def build_neighbor_map(self, data, scalar_range = None):

      number_cells = a[0].reshape(-1,3).shape[0];

      cell_ids = numpy.arange(0, number_cells, 1).reshape(-1,1)
      node_ids = a[0].reshape(-1,3);

      search_list = numpy.concatenate((cell_ids, node_ids), axis=1).astype(int)

      tmp1 = search_list[numpy.argsort(vvv[:, 1])]
      tmp2 = search_list[numpy.argsort(vvv[:, 2])]
      tmp3 = search_list[numpy.argsort(vvv[:, 3])]

      tmp44=numpy.concatenate((a1,a2), axis=1).astype(int)
      search_list=numpy.concatenate((a4,a3), axis=1).astype(int)

      for c in range(0, number_cells):
        for f in range(0, number_faces):


def connect_cells() {

    core::Matrix<size_t> fnodes = core::Matrix<size_t>::Constant(3*num_cells_, 2, 0);

    size_t j = 0;

    for (size_t i = 0; i < 3*num_cells_; i++) {

        if (i < num_cells_) {
            fnodes(i, 0) = node_ids_(j, 0);
            fnodes(i, 1) = node_ids_(j, 1);
        }
        if (i == num_cells_)
            j = 0;
        if (i >= num_cells_ && i <2*num_cells_) {
            fnodes(i, 0) = node_ids_(j, 1);
            fnodes(i, 1) = node_ids_(j, 2);
        }
        if (i == 2*num_cells_)
            j = 0;
        if (i >= 2*num_cells_ && i <3*num_cells_) {
            fnodes(i, 0) = node_ids_(j, 2);
            fnodes(i, 1) = node_ids_(j, 0);
        }
        j++;
    }

    fnodes.Sort(1);

    core::Matrix<size_t> EToE = core::Vector<size_t>::LinSpaced(num_cells_, 1, num_cells_)*core::Vector<size_t>::Ones(num_faces_per_cell_).transpose();

    core::Matrix<size_t> EToF = core::Vector<size_t>::Ones(num_cells_)*core::Vector<size_t>::LinSpaced(num_faces_per_cell_, 1, num_faces_per_cell_).transpose();

    core::Matrix<size_t> id = fnodes.col(0)*num_grid_nodes_ + fnodes.col(1) + core::Vector<size_t>::Ones(3*num_cells_);
    core::Matrix<size_t> spNodeToNode(num_faces_per_cell_*num_cells_, 4);


    spNodeToNode.col(0) = id;
    spNodeToNode.col(1) = core::Vector<size_t>::LinSpaced(num_faces_per_cell_*num_cells_, 1, num_faces_per_cell_*num_cells_).transpose();
    spNodeToNode.col(2) << core::Vector<size_t>(Eigen::Map<core::Vector<size_t>>(EToE.data(), num_faces_per_cell_*num_cells_));
    spNodeToNode.col(3) << core::Vector<size_t>(Eigen::Map<core::Vector<size_t>>(EToF.data(), num_faces_per_cell_*num_cells_));

    spNodeToNode.SortRows();

    std::vector<size_t> indices;
    for (size_t i = 0; i < num_faces_per_cell_*num_cells_-1; i++) {

        if (spNodeToNode(i, 0) == spNodeToNode(i+1, 0)) {
            indices.push_back(i);
        }
    }

    core::Matrix<size_t> matchL = core::Matrix<size_t>::Ones(indices.size()*2, 4);
    core::Matrix<size_t> matchR = core::Matrix<size_t>::Ones(indices.size()*2, 4);

    for (size_t i = 0; i < indices.size(); i++) {

        matchL.row(i) = spNodeToNode.row(indices[i]);
        matchR.row(i) = spNodeToNode.row(indices[i]+1);
    }
    size_t jj = 0;
    for (size_t i = indices.size(); i < 2*indices.size(); i++) {

        matchL.row(i) = spNodeToNode.row(indices[jj]+1);
        matchR.row(i) = spNodeToNode.row(indices[jj]);
        jj++;
    }
    for (size_t i = 0; i < 2*indices.size(); i++) {

        EToE.data()[(size_t)(matchL(i, 1)) -1] = matchR(i, 2);
        EToF.data()[(size_t)(matchL(i, 1)) -1] = matchR(i, 3);
    }

    neighbor_ids_ = core::Matrix<size_t>::Constant(num_cells_, num_faces_per_cell_, 0);
    adjacent_face_ids_ = core::Matrix<size_t>::Constant(num_cells_, num_faces_per_cell_, 0);
    core::Matrix<size_t> ones = core::Matrix<size_t>::Constant(num_cells_, num_faces_per_cell_, 1);

    neighbor_ids_ = EToE - ones;
    adjacent_face_ids_ = EToF - ones;
}


