# The MIT License (MIT)
#
# Copyright (c) 2018 Daniel Zilles
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

## @package Visualization of the grid, node points and the data.
#  This class is using VTK for visualization.
#

## Class containing all function and objects necessary for visualization:
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
class Visualization:

    ## Initialization of the class.
    #  @param self The object pointer.
    #
    def __init__(self):

        ## VTK object that provides specifications for renderers.
        ## Controls the rendering process.
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(0.8,0.8,0.8)

        ## VTK object to specify the behavior of a rendering window.
        self.render_window = vtk.vtkRenderWindow()
        self.render_window.AddRenderer(self.renderer)
        self.render_window.SetSize(1920,1080)

        ## VTK window interactor for mouse interaction
        self.render_window_interactor = vtk.vtkRenderWindowInteractor()
        self.render_window_interactor.SetRenderWindow(self.render_window)

        ##Corrdinate axes
        #self.axes = vtk.vtkAxesActor()
        #self.renderer.AddActor(self.axes);

        ## VTK object that represents a geometric structure consisting.
        ## of vertices, lines, polygons, triangle strips and the corresponding scalar/vector data. 
        self.poly_data = vtk.vtkPolyData()

        ## VTK object that maps self.poly_data to the rendering/graphics hardware/software.
        self.data_mapper = vtk.vtkPolyDataMapper()

        ## VTK object representing an object (geometry & properties) in a rendered scene. 
        ## In this case the plotted data.
        self.data_actor = vtk.vtkActor()

        ## VTK object that maps scalar values into colors via a lookup table.
        self.lookup_table = vtk.vtkLookupTable()

        ## VTK object that creates a scalar bar (legend) with labels, ticks and data values.
        self.scalar_bar_actor = vtk.vtkScalarBarActor()

        ## VTK object holding colors and their names. 
        self.named_colors = vtk.vtkNamedColors;

        ## VTK object that provides methods needed to read the data in a vtkWindow 
        ## and uses it as input to the imaging pipeline.
        #self.window_to_image_filter = vtk.vtkWindowToImageFilter()

        #self.window_to_image_filter.SetMagnification(3)
        #self.window_to_image_filter.SetInput(self.render_window)

        ## VTK object that writes PNG files.
        #self.writer = vtk.vtkPNGWriter()
        #self.writer.SetInputConnection(self.window_to_image_filter.GetOutputPort())
  
    ## Initialization of the vtk scalar bar actor.
    #  @param self The object pointer.
    #
    #  This function initializes the scalar bar actor, that renders a legend for the plotted data.
    #  Further a couloring scheme is defined in self-lookup_table.
    #
    def initialize_scalar_bar(self):

        #self.scalar_bar_actor.GetTitleTextProperty.SetFontSize(12);
        #self.scalar_bar_actor.GetTitleTextProperty.ItalicOff();
        #self.scalar_bar_actor.GetTitleTextProperty.BoldOff();
        #self.scalar_bar_actor.GetTitleTextProperty.ShadowOff();

        self.scalar_bar_actor.SetLookupTable(self.data_mapper.GetLookupTable());
        self.scalar_bar_actor.SetHeight(0.6);
        self.scalar_bar_actor.SetVerticalTitleSeparation(30);
        self.scalar_bar_actor.UnconstrainedFontSizeOn();
        self.scalar_bar_actor.SetTitle("var?");
        self.scalar_bar_actor.SetNumberOfLabels(10);
        self.lookup_table.SetTableRange(0, 1);
        self.lookup_table.SetHueRange(0, 1);
        self.lookup_table.SetSaturationRange(0, 0);
        self.lookup_table.SetValueRange(0, 1);
        self.lookup_table.Build();
        self.data_mapper.SetLookupTable(self.lookup_table);
        self.scalar_bar_actor.SetLookupTable(self.lookup_table);

        self.renderer.AddActor(self.scalar_bar_actor);

    ## Initialization of the data plot.
    #  @param self The object pointer.
    #  @param x Numpy matrix with x-coordinates, corresponding to the data points.
    #  @param y Numpy matrix with y-coordinates, corresponding to the data points.
    #
    #  This function builds a grid from the data points and should be called 
    #  every time the x and y coordinates change.
    #
    def initialize_data_plot(self, x, y):

        points = vtk.vtkPoints()
        delaunay = vtk.vtkDelaunay2D()

        cell_data = vtk.vtkDoubleArray()

        glyphFilter = vtk.vtkVertexGlyphFilter()
        point_mapper = vtk.vtkPolyDataMapper()
      
        x = x.reshape(-1)
        y = y.reshape(-1)

        for i in range(0, x.size):
            single_point = [x[i], y[i], 0]
            points.InsertNextPoint(single_point)

        self.poly_data.SetPoints(points);
        delaunay.SetInputData(self.poly_data);

        self.data_mapper.SetInputConnection(delaunay.GetOutputPort());
        self.data_actor.SetMapper(self.data_mapper);
        self.data_actor.GetProperty().EdgeVisibilityOff();

        self.renderer.AddActor(self.data_actor);

    ## Set new data and render the image.
    #  @param self The object pointer.
    #  @param data Numpy matrix with the data values.
    #  @param scalar_range Range of the data. Should be a list with 2 entries or None.
    #
    #  This function plots the data.
    #  The scalar_range is contains the range of the couloring scheme.
    #  If the scalar_range is None, the minimum and maximum number in the data matrix are chosen.
    #
    def update_data_plot(self, data, scalar_range = None):

        cell_data = vtk.vtkDoubleArray()

        cell_data = numpy_support.numpy_to_vtk(num_array=data.ravel(), deep=False, array_type=vtk.VTK_DOUBLE)

        self.poly_data.GetPointData().SetScalars(cell_data);

        #adaptiveColorBar = false;
        #double scalarRange[2] = {0.8, 1.2};

        if (scalar_range == None):
            self.data_mapper.SetScalarRange(self.poly_data.GetScalarRange())
        else:
            self.data_mapper.SetScalarRange(scalar_range)

        self.render_window.Render()

    ## Set a grid to be rendered 
    #  @param self The object pointer.
    #  @param x Numpy matrix with x-coordinates, corresponding to the data points.
    #  @param y Numpy matrix with y-coordinates, corresponding to the data points.
    #
    # A grid is build from the x- and y-coordinates. The grid will be rendered together with other
    # objects or data every time a new image is rendered.
    #
    def plot_grid(self, x, y):

        points = vtk.vtkPoints()
        poly_data = vtk.vtkPolyData()
        delaunay = vtk.vtkDelaunay2D()
        mesh_mapper = vtk.vtkPolyDataMapper()
        actor_mesh = vtk.vtkActor()
        glyph_filter = vtk.vtkVertexGlyphFilter()
        point_mapper = vtk.vtkPolyDataMapper()

        x = x.reshape(-1)
        y = y.reshape(-1)

        for i in range(0, x.size):
            single_point = [x[i], y[i], 0]
            points.InsertNextPoint(single_point)


        poly_data.SetPoints(points)
        delaunay.SetInputData(poly_data)
        mesh_mapper.SetInputConnection(delaunay.GetOutputPort())
        actor_mesh.SetMapper(mesh_mapper)
        actor_mesh.GetProperty().EdgeVisibilityOn()
        actor_mesh.GetProperty().SetColor(0, 0, 0)
        glyph_filter.SetInputData(poly_data)
        point_mapper.SetInputConnection(glyph_filter.GetOutputPort())
        point_actor = vtk.vtkActor()
        point_actor.GetProperty().SetPointSize(5)
        point_actor.SetMapper(point_mapper)
        actor_mesh.GetProperty().SetRepresentationToWireframe()
        self.renderer.AddActor(actor_mesh)

    ## Renders the current data.
    #  @param self The object pointer.
    #
    def Render(self):

        self.renderer.ResetCamera()
        self.render_window.Render()
        self.render_window_interactor.Start()

    ## Writes the current render_window to a PNG file
    #  @param self The object pointer.
    #  @param path_to_file Path to the output PNG file.
    #
    def WriteImage(self, path_to_file):
 
       self.window_to_image_filter = vtk.vtkWindowToImageFilter()
       #self.window_to_image_filter.Update()
       self.window_to_image_filter.SetInput(self.render_window)
       self.writer = vtk.vtkPNGWriter()
       self.writer.SetInputConnection(self.window_to_image_filter.GetOutputPort())
       self.writer.SetInputConnection(self.window_to_image_filter.GetOutputPort())
       self.writer.SetFileName(path_to_file);
       self.writer.Write();
