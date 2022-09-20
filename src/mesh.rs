#[allow(dead_code)]
use ndarray::prelude::*;
extern crate pyo3;
use pyo3::prelude::*;
//use numpy::ndarray::array;
use numpy::PyArray2;
use pyo3::{pymodule, types::PyModule, PyResult, Python};

#[pyclass]
/// A struct that represents a numerical mesh
pub struct Mesh
{
	/// Dimension of the space
    dimension: usize,

	/// Number of node coordinates
	number_nodes: usize,

	/// Number of cells
	number_cells: usize,

	/// Cell type (See enum CellType for details)
	cell_type: Option<CellType>,

	/// Coordinates of the grid points
	coordinates: Array2<f64>,

    // Node ids of the grid
 	node_ids : Array2<usize>,

	// Coordinates of the boundary points
	// coordinates_boundary : ArrayBase<f64, Ix3>,

	// Indices of the boundary points
	// boundary_ids

	// physical_names
}

#[pymethods]
impl Mesh {
	/// Initialize grid of variable dimension and size 
    #[new]
	pub fn new (size: usize, dim: usize) -> Mesh
	{
		let shape = (size, dim);
		Mesh
		{
			 dimension: dim,
			 number_nodes: size,
		     number_cells: 0,
		     cell_type: None,
		     coordinates: Array2::<f64>::zeros(shape),
		     node_ids: Array2::<usize>::zeros(shape),
        }
	}
		
	/// Returns the dimension
    #[getter]
	pub fn get_dimension (&self) -> usize
	{
		return self.dimension;
	}
	
	/// Returns the number of node coordinates
    #[getter]
	pub fn get_node_count (&self) -> usize
	{
		return self.number_nodes
	}
	
	/// Returns the number of cells
    #[getter]
	pub fn get_cell_count (&self) -> usize
	{
		return self.number_cells;
	}
 
	pub fn set_coordinates<'py> (&mut self, coordinates: &PyArray2<f64>, nodes: &PyArray2<usize>)
	{
        self.coordinates  = unsafe{coordinates.as_array().into_owned()};
        self.node_ids  = unsafe{nodes.as_array().into_owned()};

        self.cell_type = check_cell_type(self.node_ids.shape()[1], self.coordinates.shape()[1]).ok();
	}


	// Return reference to the coordinates
	//pub fn ref_coordinates(&self) -> &Array<DataType, Dimension> {
	//	return &self.coordinates;
	//}
    //
}

#[pymodule]
fn pazuzu(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Mesh>()?;
//    m.add_wrapped(wrap_pyfunction!(set_coordinates)).unwrap();
    Ok(())
}


enum CellType {
    Line,
	Triangle,
    Quadrangle,
    Tetrahedron,
    Hexahedron,
    Prism,
    Pyramid,
}

/// Check the cell type
fn check_cell_type(nodes: usize, dimension: usize) -> Result<CellType, String> {

    if dimension == 1 {

        match nodes {

            2 => return Ok(CellType::Line),
            _ => return Err("Cell type not implemented".to_string()),
        }
    }
    else if dimension == 2 {

        match nodes {

            3 => return Ok(CellType::Triangle),
            4 => return Ok(CellType::Quadrangle),
            _ => return Err("Cell type not implemented".to_string()),
        }
    }
    else if dimension == 3 {

        match nodes {

            4 => return Ok(CellType::Tetrahedron),
            8 => return Ok(CellType::Hexahedron),
            6 => return Ok(CellType::Prism),
            5 => return Ok(CellType::Pyramid),
            _ => return Err("Cell type not implemented".to_string()),
        }
    }
    else {
        return Err("Dimension must be 1D, 2D or 3D".to_string());
    }
}
