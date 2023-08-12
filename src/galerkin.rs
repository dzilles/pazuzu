use crate::numerical_basics;
use nalgebra::{DVector, DMatrix};
use crate::mesh::Mesh;
use libm::cos;
use libm::sin;
use std::f64::consts::PI;

/// A struct that represents a Galerkin mesh
pub struct GalerkinMesh
{
    // Order of the Lagrange polynomial
    order: usize,

    // Number of Galerkin nodes in a single cell
    //nodes_per_cell: usize,

    // Total number of Galerkin nodes
    //num_nodes: usize,

    // Number of Galerkin nodes on the face of a single cell
    //nodes_per_face: usize,

    // Coordinates for the reference nodes of the standard triangle TODO baycentric? TODO size
    //r_1: DVector<f64>, 
    //r_2: DVector<f64>, 
    //r_3: DVector<f64>,

    // Vandermonde matrix of size TODO
    //vandermonde: DMatrix<f64>,

    // Inverse of vandermonde matrix of size [m_noNodes, m_noNodes]
    //inverse_vandermonde: DMatrix<f64>,

    // Mass matrix of size [m_noNodes, m_noNodes]
    //mass_matrix, DMatrix<f64>;

    // Differentiation matrices of size [m_noNodes, m_noNodes]
    //diff_r_1: DMatrix<f64>, 
    //diff_r_2: DMatrix<f64>,
    //diff_r_3: DMatrix<f64>,

    // Global coordinates of all nodes, matrix of size [m_noNodes, m_noCells] 
    //x_1: DMatrix<f64>,
    //x_2: DMatrix<f64>,
    //x_3: DMatrix<f64>,

    // indices to the nodes at the cell boundaries, matrix of size [m_noFaces*m_nodesFaces, m_noCells]
    //face_node_ids: DMatrix<usize>;
    // TODO Remove m_fmask, ..., m_fmask1
    //m_fmask: DVector<usize>

    // Node coordinates on the faces, matrix of size [m_noFaces*m_nodesFace, m_noCells]
    //x_1_face: DMatrix<f64>, 
    //x_2_face: DMatrix<f64>, 
    //x_3_face: DMatrix<f64>,

    // Surface integral, matrix of size [m_noNodes, m_noFaces*m_noNodesFace]
    //surface_integral:DMatrix<f64>,

    // Metric constants, matrix of size [m_noNodes, m_noCells]
    //rx_11: DMatrix<f64>, 
    //rx_12: DMatrix<f64>, 
    //rx_13: DMatrix<f64>, 
    //rx_21: DMatrix<f64>, 
    //rx_22: DMatrix<f64>, 
    //rx_23: DMatrix<f64>, 
    //rx_31: DMatrix<f64>, 
    //rx_32: DMatrix<f64>, 
    //rx_33: DMatrix<f64>, 

    // Volume Jacobian, matrix of size [m_noNodes, m_noCells]
    //jacobian: DMatrix<f64>,

    // Outward normals on the nodes on the cell faces, matrix of size [m_noFaces*m_noNodesFace, m_noCells]
    //normals_x_1: DMatrix<f64>, 
    //normals_x_2: DMatrix<f64>, 
    //normals_x_3: DMatrix<f64>, 

    // Volume Jacobian, matrix of size [m_noFaces*m_noNodesFace, m_noCells]
    //surface_jacobian: DMatrix<f64>,

    // TODO comment
    //fscale: DMatrix<f64>,

    // TODO comment
    //weak_diff_r_1: DMatrix<f64>,
    //weak_diff_r_2: DMatrix<f64>,
    //weak_diff_r_3: DMatrix<f64>,

    // Coordinates of the grid vertices, vector of size [m_noVertices]
    //core::Vector<p_float> grid_x_, grid_y_, grid_z_;

    // Vertex identifier vertexId(cid,vid) of vertex vid in cell cid, size [m_noCells, m_noFaces]
    //node_ids: DMatrix<usize>,

    // Neighbor cell identifier neighborId(cid,fid) of face fid in cell cid, size [m_noCells, m_noFaces]
    //neighbor_ids: DMatrix<usize>,

    // Neighbor local face identifier neighborId(cid,fid) of face fid in cell cid, size [m_noCells, m_noFaces]
    //adjacent_face_ids: DMatrix<usize>,

    // Boundary Condition Id for ...
    //boundary_condition_ids: DMatrix<usize>,
}

impl GalerkinMesh {

    pub fn new(order:usize) -> GalerkinMesh {

        let tmp = GalerkinMesh {
            order: order,
        };

        Self::alpha_optimized_nodes(order);

        let b = numerical_basics::jacobi_gauss_lobatto(0.0, 0.0, 5);
        tmp
    }

    //pub fn new(msh: Mesh, order:usize) {


        //(x_1_hat, x_2_hat, x_3_hat) = alpha_optimized_nodes(order);

        //(r_1, r_2, r_3) = standard_triangle(x_1_hat, x_2_hat, x_3_hat);

        //vandermonde = vandermonde(order, r_1, r_2, r_3);
        //inverse_vandermonde = vandermonde.inverse();
        //mass_matrix = inverse_vandermonde.transpose()*inverse_vandermonde;

    //}
    fn alpha_optimized_nodes(order: usize) -> (DVector<f64>, DVector<f64>, DVector<f64>){ 

        let number_nodes = (order+1)*(order+2)/2;

        let mut tmp_x: DVector<f64> = DVector::from_element(number_nodes, 0.0);
        let mut tmp_y: DVector<f64> = DVector::from_element(number_nodes, 0.0);
        let mut tmp_z: DVector<f64> = DVector::from_element(number_nodes, 0.0);

        let mut alpha = 5.0/3.0;
        let alpha_table = [
            0.0000, 0.0000, 1.4152, 0.1001, 0.2751,
            0.9800, 1.0999, 1.2832, 1.3648, 1.4773,
            1.4959, 1.5743, 1.5770, 1.6223, 1.6258];

        if order < 15 {
            alpha = alpha_table[order-1];
        }

        let mut l1: DVector<f64> = DVector::from_element(number_nodes, 0.0);
        let mut l2: DVector<f64> = DVector::from_element(number_nodes, 0.0);
        let mut l3: DVector<f64> = DVector::from_element(number_nodes, 0.0);

        let mut cnt: usize = 0;
        for n in 1..order+2 {
            for m in 1..order+3-n {

                l1[cnt] = (n as f64 - 1.0)/order as f64;
                l3[cnt] = (m as f64 - 1.0)/order as f64;
                l2[cnt] = 1.0 - l1[cnt] - l3[cnt];
                cnt+=1;
            }
        }

        tmp_x = &l3 - &l2;
        tmp_y = (2.0*&l1 - &l2 - &l3)/f64::sqrt(3.0);

        let mut blend1: DVector<f64> = 4.0 * l2.component_mul(&l3);
        let mut blend2: DVector<f64> = 4.0 * l1.component_mul(&l3);
        let mut blend3: DVector<f64> = 4.0 * l1.component_mul(&l2);

        let mut lw1: DVector<f64> = &l3 - &l2;
        let mut lw2: DVector<f64> = &l1 - &l3;
        let mut lw3: DVector<f64> = &l2 - &l1;

        let mut warpf1: DVector<f64> = numerical_basics::warpfactor(order, &lw1);
        let mut warpf2: DVector<f64> = numerical_basics::warpfactor(order, &lw2);
        let mut warpf3: DVector<f64> = numerical_basics::warpfactor(order, &lw3);


  
        let mut warp1 = blend1.component_mul(&warpf1).component_mul( &(alpha*alpha*l1.component_mul(&l1)).add_scalar(1.0) );


        let mut warp2 = blend2.component_mul(&warpf2).component_mul( &(alpha*alpha*l2.component_mul(&l2)).add_scalar(1.0) );
        let mut warp3 = blend3.component_mul(&warpf3).component_mul( &(alpha*alpha*l3.component_mul(&l3)).add_scalar(1.0) );

        tmp_x = tmp_x + 1.0*(&warp1) + cos(2.0*PI/3.0)*(&warp2) + cos(4.0*PI/3.0)*(&warp3);
        tmp_y = tmp_y + 0.0*warp1 + sin(2.0*PI/3.0)*warp2 + sin(4.0*PI/3.0)*warp3;

        (tmp_x, tmp_y, tmp_z)
    }


    fn standard_triangle(x: &DVector<f64>, y: &DVector<f64>) -> (DVector<f64>, DVector<f64>) {

        let l1 = ((f64::sqrt(3.0)*y).add_scalar(1.0))/3.0;
        let l2 = (-3.0*x - (f64::sqrt(3.0)*y).add_scalar(2.0))/6.0;
        let l3 = (3.0*x - (f64::sqrt(3.0)*y).add_scalar(2.0))/6.0;

        let r = &l3 - &l1 - &l2;
        let s = l1 - l2 - l3;

        (r, s)
    }

    fn simplex(a: &DVector<f64>, b: &DVector<f64>, order_1: usize, order_2: usize) -> DVector<f64>{

        let h1 = numerical_basics::jacobi_polynomial(0.0, 0.0, order_1, &a);
        let h2 = numerical_basics::jacobi_polynomial(2.0*(order_1 as f64)+1.0, 0.0, order_2, &b);

        let mut h3 = (-b).add_scalar(1.0);

        for x in h3.iter_mut() {

            *x = x.powf(order_1 as f64);
        }

        let p = f64::sqrt(2.0)*h1.component_mul(&h2).component_mul(&h3);

        p
    }


    fn vandermonde(order: usize, r: &DVector<f64>, s: &DVector<f64>) {

        let size = r.shape().0;

        let vandermonde: DMatrix<f64> = DMatrix::from_element(size, (order+1)*(order+2)/2, 0.0);
        let (a, b) = numerical_basics::orthonormalize(r, s);
    }


    //fn diff_matrices() {
    //}

    //fn warpfactor(usize order, DVector<f64> r_out) {


    //}

/*
    }

    fn vandermonde_1d() {
    }

    fn grad_vandermonde() {
    }

    fn grad_simplex() {
    }

    )*fn surf_integral() {
    }*/
}

#[cfg(test)]
mod tests {
    use crate::galerkin::*;
    use std::fs::{File, read_to_string};
    use std::io::{BufWriter, Write};
    use serde_json::json;
    use serde::{Deserialize, Serialize};

    #[cfg(feature = "generate_test_data")]
    fn generate_test_data_alpha_optimized_nodes() {

        let (r, s, t) = GalerkinMesh::alpha_optimized_nodes(8);

        #[derive(Serialize)]
        struct TestData {
            r: Vec<f64>,
            s: Vec<f64>,
            t: Vec<f64>,
        }

        let mut r: Vec<f64> = r.data.into();
        let mut s: Vec<f64> = s.data.into();
        let mut t: Vec<f64> = t.data.into();

        let file = File::create("tests/data/unit_tests/test_alpha_optimized_nodes.json");
        let mut writer = BufWriter::new(file.expect("Failed to create file"));
        let test_data = json!({
            "r": r,
            "s": s,
            "t": t,
        });
        serde_json::to_writer(&mut writer, &test_data);

        //TODO close file
    }

    #[test] 
    fn alpha_optimized_nodes() {

        #[cfg(feature = "generate_test_data")]
        generate_test_data_alpha_optimized_nodes();

        let (r, s, t) = GalerkinMesh::alpha_optimized_nodes(8);

        #[derive(Deserialize)]
        struct TestData {
            r: Vec<f64>,
            s: Vec<f64>,
            t: Vec<f64>,
        }

        let mut r: Vec<f64> = r.data.into();
        let mut s: Vec<f64> = s.data.into();
        let mut t: Vec<f64> = t.data.into();

        let file_content = read_to_string("tests/data/unit_tests/test_alpha_optimized_nodes.json").expect("Error reading file");

        let test_data: TestData = serde_json::from_str(&file_content).expect("Error creating test data");

        let eps = 10E-12;
     
        assert!(r.len() == s.len());
        assert!(s.len() == t.len());

        for i in 0..r.len() {
           
            float_eq::assert_float_eq!(r[i], test_data.r[i], abs <= eps);
            float_eq::assert_float_eq!(s[i], test_data.s[i], abs <= eps);
            float_eq::assert_float_eq!(t[i], test_data.t[i], abs <= eps);
        }
    }

    #[test] 
    fn standard_triangle() {

        assert!(false);
    }

    #[test] 
    fn simplex() {

        assert!(false);
    }

    #[test] 
    fn vandermonde() {

        assert!(false);
    }

    #[test] 
    fn orthonormalize() {

        assert!(false);
    }

    #[test] 
    fn diff_matrices() {

        assert!(false);
    }

    #[test] 
    fn global_node_coorinates() {

        assert!(false);
    }

    #[test] 
    fn surface_integral() {

        assert!(false);
    }
    #[test] 
    fn geometric_factors() {

        assert!(false);
    }

    #[test] 
    fn normals() {

        assert!(false);
    }

    #[test] 
    fn weak_terms() {

        assert!(false);
    }
}
