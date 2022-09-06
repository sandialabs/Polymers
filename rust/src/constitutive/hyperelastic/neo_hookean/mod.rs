use ndarray::Array2;
use crate::constitutive::IDENTITY_TENSOR;
use crate::constitutive::get_left_cauchy_green_tensor;
use crate::constitutive::hyperelastic::Hyperelastic;

pub mod test;

pub struct NeoHookean {
    pub hyperelastic: Hyperelastic
}

impl NeoHookean {
    pub fn init() -> NeoHookean {
        NeoHookean {
            hyperelastic: Hyperelastic::init()
        }
    }
    pub fn get_cauchy_stress(&self, deformation_gradient: &Array2<f64>) -> Array2<f64> {
        let left_cauchy_green_tensor =& get_left_cauchy_green_tensor(deformation_gradient);
        left_cauchy_green_tensor - &*IDENTITY_TENSOR
    }
}
