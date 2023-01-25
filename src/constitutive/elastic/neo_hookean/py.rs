use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let neo_hookean = PyModule::new(py, "neo_hookean")?;
    parent_module.add_submodule(neo_hookean)?;
    neo_hookean.add_class::<NeoHookean>()?;
    Ok(())
}

/// The Neo-Hookean model.
#[pyclass]
#[derive(Copy, Clone)]
pub struct NeoHookean
{
    /// The shear modulus in units of ???.
    #[pyo3(get)]
    pub shear_modulus: f64,

    /// The bulk modulus in units of ???.
    #[pyo3(get)]
    pub bulk_modulus: f64
}

#[pymethods]
impl NeoHookean
{
    #[new]
    pub fn init(shear_modulus: f64, bulk_modulus: f64) -> Self
    {
        NeoHookean
        {
            shear_modulus,
            bulk_modulus
        }
    }
    /// The true stress as a function of the deformation gradient.
    pub fn true_stress(&self, deformation_gradient: [[f64; 3]; 3]) -> PyResult<[[f64; 3]; 3]>
    {
        Ok(super::NeoHookean::init(self.shear_modulus, self.bulk_modulus).true_stress(&deformation_gradient))
    }
}