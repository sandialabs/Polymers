use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let fjc = PyModule::new(py, "fjc")?;
    super::thermodynamics::py::register_module(py, fjc)?;
    parent_module.add_submodule(fjc)?;
    fjc.add_class::<FJC>()?;
    Ok(())
}

/// The freely-jointed chain (FJC) model.
#[pyclass]
pub struct FJC
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
    pub thermodynamics: super::thermodynamics::py::FJC
}

#[pymethods]
impl FJC
{
    #[new]
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            thermodynamics: super::thermodynamics::py::FJC::init(number_of_links, link_length, hinge_mass)
        }
    }
}
