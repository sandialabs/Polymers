use pyo3::prelude::*;
use numpy::
{
    IntoPyArray,
    PyArrayDyn,
    PyReadonlyArrayDyn
};
use crate::physics::BOLTZMANN_CONSTANT;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let isotensional = PyModule::new(py, "isotensional")?;
    super::asymptotic::py::register_module(py, isotensional)?;
    super::legendre::py::register_module(py, isotensional)?;
    parent_module.add_submodule(isotensional)?;
    isotensional.add_class::<MORSEFJC>()?;
    Ok(())
}

/// The Morse link potential freely-jointed chain (Morse-FJC) model thermodynamics in the isotensional ensemble.
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

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using an asymptotic approach.
    #[pyo3(get)]
    pub asymptotic: super::asymptotic::py::MORSEFJC,

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using a Legendre transformation.
    #[pyo3(get)]
    pub legendre: super::legendre::py::MORSEFJC
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
            asymptotic: super::asymptotic::py::MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy),
            legendre: super::legendre::py::MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        }
    }
}