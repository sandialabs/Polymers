use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let thermodynamics = PyModule::new(py, "thermodynamics")?;
    super::isometric::py::register_module(py, thermodynamics)?;
    super::isotensional::py::register_module(py, thermodynamics)?;
    parent_module.add_submodule(thermodynamics)?;
    thermodynamics.add_class::<MORSEFJC>()?;
    Ok(())
}

/// The Morse link potential freely-jointed chain (Morse-FJC) model thermodynamics.
#[pyclass]
#[derive(Copy, Clone)]
pub struct MORSEFJC
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

    /// The energy of each link in the chain in units of J/mol.
    #[pyo3(get)]
    pub link_energy: f64,
    
    /// The thermodynamic functions of the model in the isometric ensemble.
    #[pyo3(get)]
    pub isometric: super::isometric::py::MORSEFJC,
    
    /// The thermodynamic functions of the model in the isotensional ensemble.
    #[pyo3(get)]
    pub isotensional: super::isotensional::py::MORSEFJC
}

#[pymethods]
impl MORSEFJC
{
    #[new]
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64, link_energy: f64) -> Self
    {
        MORSEFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            link_energy,
            isometric: super::isometric::py::MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy),
            isotensional: super::isotensional::py::MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        }
    }
}
