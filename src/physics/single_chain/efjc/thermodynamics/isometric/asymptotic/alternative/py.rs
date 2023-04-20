use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let alternative = PyModule::new(py, "alternative")?;
    super::legendre::py::register_module(py, alternative)?;
    parent_module.add_submodule(alternative)?;
    alternative.add_class::<EFJC>()?;
    Ok(())
}

/// The extensible freely-jointed chain (EFJC) model thermodynamics in the isometric ensemble approximated using a alternative asymptotic approach.
#[pyclass]
#[derive(Copy, Clone)]
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

    /// The thermodynamic functions of the model in the isometric ensemble approximated using a alternative asymptotic approach and a Legendre transformation.
    #[pyo3(get)]
    pub legendre: super::legendre::py::EFJC
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
            legendre: super::legendre::py::EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness)
        }
    }
}