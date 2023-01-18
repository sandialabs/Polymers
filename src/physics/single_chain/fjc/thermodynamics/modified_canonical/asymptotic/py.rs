use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let asymptotic = PyModule::new(py, "asymptotic")?;
    super::weak_potential::py::register_module(py, asymptotic)?;
    super::strong_potential::py::register_module(py, asymptotic)?;
    parent_module.add_submodule(asymptotic)?;
    asymptotic.add_class::<FJC>()?;
    Ok(())
}

/// The freely-jointed chain (FJC) model thermodynamics in the modified canonical ensemble approximated using an asymptotic approach.
#[pyclass]
#[derive(Copy, Clone)]
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

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using an asymptotic approach valid for weak potentials.
    #[pyo3(get)]
    pub weak_potential: super::weak_potential::py::FJC,

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using an asymptotic approach valid for strong potentials.
    #[pyo3(get)]
    pub strong_potential: super::strong_potential::py::FJC
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
            weak_potential: super::weak_potential::py::FJC::init(number_of_links, link_length, hinge_mass),
            strong_potential: super::strong_potential::py::FJC::init(number_of_links, link_length, hinge_mass)
        }
    }
}