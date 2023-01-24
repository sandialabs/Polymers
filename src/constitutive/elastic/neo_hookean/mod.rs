#[cfg(feature = "python")]
pub mod py;

use crate::constitutive::
{
    jacobian,
    deviatoric_incompressible_left_cauchy_green
};

/// The structure of the Neo-Hookean model.
pub struct NeoHookean
{
    /// The shear modulus in units of ???.
    pub shear_modulus: f64,

    /// The bulk modulus in units of ???.
    pub bulk_modulus: f64
}

/// The implemented functionality of the Neo-Hookean model.
impl NeoHookean
{
    /// Initializes and returns an instance of the Neo-Hookean model.
    pub fn init(shear_modulus: f64, bulk_modulus: f64) -> Self
    {
        NeoHookean
        {
            shear_modulus,
            bulk_modulus
        }
    }
    /// The true stress as a function of the deformation gradient.
    pub fn true_stress(&self, deformation_gradient: &[[f64; 3]; 3]) -> [[f64; 3]; 3]
    {
        let jac = jacobian(deformation_gradient);
        let dev_b_bar = deviatoric_incompressible_left_cauchy_green(deformation_gradient);
        let fac = self.shear_modulus/jac;
        let sph = self.bulk_modulus*(jac - 1.0);
        [
            [fac*dev_b_bar[0][0] + sph, fac*dev_b_bar[0][1], fac*dev_b_bar[0][2]],
            [fac*dev_b_bar[1][0], fac*dev_b_bar[1][1] + sph, fac*dev_b_bar[1][2]],
            [fac*dev_b_bar[2][0], fac*dev_b_bar[2][1], fac*dev_b_bar[2][2] + sph]
        ]
    }
}