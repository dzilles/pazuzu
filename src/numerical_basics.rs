// The MIT License (MIT)
//
// Copyright (c) 2018 Daniel Zilles
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

use nalgebra::{DVector, DMatrix};
use gauss_quad::GaussJacobi;
use float_eq::assert_float_eq;
use libm::tgamma;

/// TODO
pub fn jacobi_gauss_lobatto(alpha: f64, beta: f64, order: usize) -> DVector<f64> {

    let size = order + 1;
    let mut return_value: DVector<f64> = DVector::from_element(size, 0.0);

    if order == 1 {
        
        return_value[0] = -1.0;
        return_value[1] = 1.0;
    } 
    else {
 
        let (nodes, weights) = GaussJacobi::nodes_and_weights(order-1, alpha + 1.0, beta + 1.0);
        let mut sorted_nodes: DVector<f64> = DVector::from_element(order-1, 0.0);
        return_value[0] = -1.0;
        return_value[order] = 1.0;

        let loopcnt = (order-1)/2;

        for i in 0..loopcnt {

            let mut shift = 0;

            if ((order-1)%2 != 0) && ((i as f64)/((loopcnt -1) as f64) > 0.5) {
                shift = 1
            }

            sorted_nodes[2*i+shift] = nodes[order-2-i];
            sorted_nodes[2*i+1+shift] = nodes[i];
        }
        for i in 1..order {

            return_value[i] = sorted_nodes[i-1];
        }

    }
    return_value
}

pub fn vandermonde_1d(order: usize, x: &DVector<f64>) -> DMatrix<f64>{


    let size = x.shape().0;
    let mut vandermonde_1d = DMatrix::from_element(size, order + 1, 0.0);

    for col_idx in 0..order+1 {

        let jacobi = jacobi_polynomial(0.0, 0.0, col_idx, &x);
        for row_idx in 0..size {


            vandermonde_1d[(row_idx, col_idx)] = jacobi[row_idx];
        }
    }
    vandermonde_1d
}

pub fn jacobi_polynomial(alpha: f64, beta: f64, order: usize, x: &DVector<f64>) -> DVector<f64>{

    let size = x.shape().0;
    let mut polynomial_components = DMatrix::from_element(order+1, size, 0.0);
    let tmp = alpha + beta + 1.0;
    let gamma0 = libm::pow(2.0, tmp)/(tmp)*libm::tgamma(alpha + 1.0)*libm::tgamma(beta + 1.0)/libm::tgamma(tmp);
    
    polynomial_components.row_mut(0).add_scalar_mut(1.0/libm::sqrt(gamma0));

    if order > 0 {
    
        let gamma1 = (alpha + 1.0)*(beta + 1.0)/(alpha + beta + 3.0)*gamma0;
        polynomial_components.row_mut(1).add_scalar_mut((alpha-beta)/2.0);


        for (mut row, coeff) in polynomial_components.row_mut(1).iter_mut().zip(x.iter()) {
            *row += *coeff*(alpha + beta +2.0)/2.0;
            *row /= libm::sqrt(gamma1);
        }
    }

    if order > 1 {

        let mut aold: f64 = 2.0/(2.0+alpha+beta)*libm::sqrt((alpha+1.0)*(beta+1.0)/(alpha+beta+3.0));

        for row_idx in 1..order {

            let i = row_idx as f64;
            let h1 = 2.0*i+alpha+beta;
            let anew = 2.0/(h1+2.0)*libm::sqrt((i+1.0)*(i+1.0+alpha+beta)*(i+1.0+alpha)*(i+1.0+beta)/(h1+1.0)/(h1+3.0));
            let bnew = - (alpha*alpha-beta*beta)/h1/(h1+2.0);
 
            for col_idx in 0..size {
                polynomial_components[(row_idx+1, col_idx)] = 1.0/anew*(-aold*polynomial_components[(row_idx-1,col_idx)] + (x[col_idx]-bnew)*polynomial_components[(row_idx, col_idx)]);

            }
            aold = anew;
        }
    }
    DVector::from_iterator(size, polynomial_components.row_mut(order).iter().cloned())
}

pub fn legendre_polynomial(order: usize, x: &DVector<f64>) -> DVector<f64>{

    jacobi_polynomial(0.0, 0.0, order, x)
}

pub fn warpfactor(order: usize, x: &DVector<f64>) -> DVector<f64>{

    let lglr = jacobi_gauss_lobatto(0.0, 0.0, order);

    let half_order: f64 = (order as f64)/2.0;

    let req: DVector<f64> = DVector::from_iterator(order+1, (0..order+1).map(|num| (num as f64-half_order)/half_order));

    let veq = vandermonde_1d(order, &req);

    let nr: usize = x.shape().0;

    let mut pmat: DMatrix<f64> = DMatrix::from_element(order+1, nr, 0.0);

    for row_idx in 0..order+1 {

        let jacobi = jacobi_polynomial(0.0, 0.0, row_idx, &x);

        for col_idx in 0..nr {

            pmat[(row_idx, col_idx)] = jacobi[col_idx];

        }
    }

    let decomposition = veq.transpose().lu();

    let lmat = decomposition.solve(&pmat).expect("Linear resolution failed");

    let mut warp: DVector<f64> = DVector::from_element(order+1, 0.0);
    warp = lmat.transpose()*(lglr.clone() - req.clone());

    let mut zeroof: DVector<f64> = DVector::from_element(x.shape().0, 0.0);
    for i in 0..x.shape().0 {

        if x.abs()[i] < 1.0
        {
            zeroof[i] = 1.0
        }
        else
        {
            zeroof[i] = 0.0
        }
    }

    let mut sf = zeroof.component_mul(&x);
    let mut sf = sf.component_mul(&sf);
    let mut sf = DVector::from_element(x.shape().0, 1.0) - sf;

    zeroof.add_scalar_mut(-1.0);

    let warp1 = warp.component_div(&sf);
    let warp2 = warp.component_mul(&zeroof);

    warp1+warp2
}

pub fn orthonormalize(r: &DVector<f64>, s: &DVector<f64>) -> (DVector<f64>, DVector<f64>) {

    let mut a = DVector::from_element(r.shape().0, 0.0);

    for n in 0..r.shape().0 {

        if s[n] != 1.0
        {
            a[n] = 2.0*(1.0+r[n])/(1.0-s[n])-1.0;
        }
        else
        {
            a[n] = -1.0;
        }
    }
    (a,s.clone())
}

#[cfg(test)]
mod tests {
    use crate::numerical_basics;
    use nalgebra::{DVector, DMatrix};

    #[test] 
    fn jacobi_gauss_lobatto() {

        let nodes = numerical_basics::jacobi_gauss_lobatto(0.0, 0.0, 5);
        let eps = 10E-12;
        float_eq::assert_float_eq!(nodes[0], -1.0, abs <= eps);
        float_eq::assert_float_eq!(nodes[1], 0.7650553239294648, abs <= eps);
        float_eq::assert_float_eq!(nodes[2], -0.7650553239294647, abs <= eps);
        float_eq::assert_float_eq!(nodes[3], 0.2852315164806454, abs <= eps);
        float_eq::assert_float_eq!(nodes[4], -0.2852315164806451, abs <= eps);
        float_eq::assert_float_eq!(nodes[5], 1.0, abs <= eps);
    }

    #[test] 
    fn jacobi_polynomial() {

        let vec: DVector<f64> = DVector::from_iterator(4, [0.0, 1.0, 2.0, 3.0].iter().cloned());
        let eps = 10E-12;

        let jacobi0 = numerical_basics::jacobi_polynomial(1.0, 2.0, 0, &vec);

        let jacobi = numerical_basics::jacobi_polynomial(1.0, 2.0, 5, &vec);
        float_eq::assert_float_eq!(jacobi[0], -0.5063078670631143, abs <= eps);
        float_eq::assert_float_eq!(jacobi[1], 6.4807406984078595, abs <= eps);
        float_eq::assert_float_eq!(jacobi[2], 900.5191723584549, abs <= eps);
        float_eq::assert_float_eq!(jacobi[3], 8395.799574787385, abs <= eps);
    }

    #[test] 
    fn vandermonde_1d() {

        let vec: DVector<f64> = DVector::from_iterator(4, [0.0, 1.0, 2.0, 3.0].iter().cloned());
        let vandermonde = numerical_basics::vandermonde_1d(5, &vec);
        let eps = 10E-12;
        
        float_eq::assert_float_eq!(vandermonde[(0,0)], 0.7071067811865475, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(0,1)], 0.0, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(0,2)], -0.7905694150420949, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(0,3)], 0.0, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(0,4)], 0.7954951288348658, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(0,5)], 0.0, abs <= eps);
        
        float_eq::assert_float_eq!(vandermonde[(1,0)], 0.7071067811865475, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(1,1)], 1.224744871391589, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(1,2)], 1.58113883008419, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(1,3)], 1.8708286933869718, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(1,4)], 2.121320343559644, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(1,5)], 2.3452078799117166, abs <= eps);

        float_eq::assert_float_eq!(vandermonde[(2,0)], 0.7071067811865475, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(2,1)], 2.449489742783178, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(2,2)], 8.696263565463044, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(2,3)], 31.804087787578506, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(2,4)], 117.46811402461522, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(2,5)], 435.6223636936011, abs <= eps);

        float_eq::assert_float_eq!(vandermonde[(3,0)], 0.7071067811865475, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(3,1)], 3.6742346141747673, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(3,2)], 20.55480479109447, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(3,3)], 117.86220768337918, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(3,4)], 680.9438302826453, abs <= eps);
        float_eq::assert_float_eq!(vandermonde[(3,5)], 3946.9848618914166, abs <= eps);
    }

    #[test] 
    fn warpfactor() {

        let vec: DVector<f64> = DVector::from_iterator(4, [0.0, 1.0, 2.0, 3.0].iter().cloned());
        let warpf = numerical_basics::warpfactor(5, &vec);
        let eps = 10E-12;
        float_eq::assert_float_eq!(warpf[0], -0.3280045558732164 , abs <= eps);
        float_eq::assert_float_eq!(warpf[1], 0.0, abs <= eps);
        float_eq::assert_float_eq!(warpf[2], 0.0, abs <= eps);
        float_eq::assert_float_eq!(warpf[3], 0.0, abs <= eps);
    }

    #[test] 
    fn orthonormalize() {

        assert!(false)
    }
}
