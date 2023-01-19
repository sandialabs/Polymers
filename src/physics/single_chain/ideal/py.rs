use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let ideal = PyModule::new(py, "ideal")?;
    super::thermodynamics::py::register_module(py, ideal)?;
    parent_module.add_submodule(ideal)?;
    ideal.add_class::<Ideal>()?;
    Ok(())
}

/// The ideal chain model.
#[pyclass]
pub struct Ideal
{
    /// The mass of each hinge in the chain in units of kg/mol.
    #[pyo3(get)]
    pub hinge_mass: f64,
    
    /// The length of each link in the chain in units of nm.
    #[pyo3(get)]
    pub link_length: f64,
    
    /// The number of links in the chain.
    #[pyo3(get)]
    pub number_of_links: u8,
    
    /// The thermodynamic functions of the model.
    #[pyo3(get)]
    pub thermodynamics: super::thermodynamics::py::Ideal
}

#[pymethods]
impl Ideal
{
    #[new]
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        Ideal
        {
            hinge_mass,
            link_length,
            number_of_links,
            thermodynamics: super::thermodynamics::py::Ideal::init(number_of_links, link_length, hinge_mass)
        }
    }
}
