use ndarray::{
    Array2,
    linalg::{
        general_mat_mul
    }
};
use lazy_static::lazy_static;

pub mod test;
pub mod hyperelastic;

lazy_static! {
    static ref ZERO_TENSOR: Array2::<f64> = Array2::<f64>::zeros((3, 3));
    static ref IDENTITY_TENSOR: Array2::<f64> = Array2::<f64>::eye(3);
}

pub struct Constitutive {
}

impl Constitutive {
    pub fn init() -> Constitutive {
        Constitutive {
        }
    }
}

pub fn matrix_multiply(a: &Array2<f64>, b: &Array2<f64>) -> Array2<f64> {
    let mut c = Array2::<f64>::zeros((3, 3));
    general_mat_mul(1.0, a, b, 0.0, &mut c);
    c
}

pub fn get_left_cauchy_green_tensor(deformation_gradient: &Array2<f64>) -> Array2<f64> {
    matrix_multiply(&deformation_gradient.to_owned().reversed_axes(), deformation_gradient)
}
