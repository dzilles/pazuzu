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

	/// Number of nodes. Nodes are the coordinate points of a mesh
	number_nodes: usize,

	/// Number of cells. A cell is a space that is bounded by a certain number of nodes.
	number_cells: usize,

	//number_vertices: usize,

	//number_vertices_boundary: usize,

	//number_faces: usize,

	/// Cell type (See enum CellType for details). The cell type determines the number of nodes of
    /// a single cell
	cell_type: Option<CellType>,

	/// coordinates(nid, dim) x (dim=0), y (dim=1), z (dim=2) coordinates of node nid
	coordinates: Array2<f64>,

    // node_ids(vid, cid) of vertex vid in cell cid, size [number_cells, number_faces]
 	node_ids : Array2<usize>,

	// Indices of the boundary points
 	//boundary_ids : Array2<usize>,

 	//neighbor_ids : Array2<usize>,
	// physical_names
 	//adjacent_face_ids : Array2<usize>,
}

#[pymethods]
impl Mesh {
    #[new]
	/// Initialize grid of variable dimension and size 
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
	pub fn get_dimension (&self) -> usize
	{
		return self.dimension;
	}
	
	/// Returns the number of node coordinates
	pub fn get_node_count (&self) -> usize
	{
		return self.number_nodes
	}
	
	/// Returns the number of cells
	pub fn get_cell_count (&self) -> usize
	{
		return self.number_cells;
	}
 
	pub fn setup_mesh<'py> (&mut self, coordinates: &PyArray2<f64>, nodes: &PyArray2<usize>)
	{
        self.coordinates  = unsafe{coordinates.as_array().into_owned()};
        self.node_ids  = unsafe{nodes.as_array().into_owned()};

        self.cell_type = check_cell_type(self.node_ids.shape()[1], self.coordinates.shape()[1]).ok();
	}

	pub fn setup_boundary (&mut self, _boundary_nodes: &PyArray2<usize>, _group_name: String)
	{
	}
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
    // Prism,
    // Pyramid,
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
            // 6 => return Ok(CellType::Prism),
            // 5 => return Ok(CellType::Pyramid),
            _ => return Err("Cell type not implemented".to_string()),
        }
    }
    else {
        return Err("Dimension must be 1D, 2D or 3D".to_string());
    }
}

#[test]
fn test_check_cell_type() {

    assert!(check_cell_type(3, 0).is_err());
    assert!(check_cell_type(3, 4).is_err());
    assert!(check_cell_type(3, 1).is_err());
    assert!(check_cell_type(5, 2).is_err());
    assert!(check_cell_type(1, 3).is_err());
    assert!(matches!(check_cell_type(2, 1).unwrap(), CellType::Line));
    assert!(matches!(check_cell_type(3, 2).unwrap(), CellType::Triangle));
    assert!(matches!(check_cell_type(4, 2).unwrap(), CellType::Quadrangle));
    assert!(matches!(check_cell_type(4, 3).unwrap(), CellType::Tetrahedron));
    assert!(matches!(check_cell_type(8, 3).unwrap(), CellType::Hexahedron));
}
