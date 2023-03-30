use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let log_squared = PyModule::new(py, "log_squared")?;
    super::thermodynamics::py::register_module(py, log_squared)?;
    parent_module.add_submodule(log_squared)?;
    log_squared.add_class::<LOGSQUAREDFJC>()?;
    Ok(())
}

/// The log-squared link potential freely-jointed chain (log-squared-FJC) model.
#[pyclass]
pub struct LOGSQUAREDFJC
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

    /// The stiffness of each link in the chain in units of J/(mol⋅nm^2).
    #[pyo3(get)]
    pub link_stiffness: f64,
    
    /// The thermodynamic functions of the model.
    #[pyo3(get)]
    pub thermodynamics: super::thermodynamics::py::LOGSQUAREDFJC
}

#[pymethods]
impl LOGSQUAREDFJC
{
    #[new]
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        LOGSQUAREDFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            thermodynamics: super::thermodynamics::py::LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness)
        }
    }
}
