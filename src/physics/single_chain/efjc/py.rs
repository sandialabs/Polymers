use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let efjc = PyModule::new(py, "efjc")?;
    super::thermodynamics::py::register_module(py, efjc)?;
    parent_module.add_submodule(efjc)?;
    efjc.add_class::<EFJC>()?;
    Ok(())
}

/// The extensible freely-jointed chain (EFJC) model.
#[pyclass]
pub struct EFJC
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
    pub thermodynamics: super::thermodynamics::py::EFJC
}

#[pymethods]
impl EFJC
{
    #[new]
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        EFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            thermodynamics: super::thermodynamics::py::EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness)
        }
    }
}
