use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let thermodynamics = PyModule::new(py, "thermodynamics")?;
    super::isotensional::py::register_module(py, thermodynamics)?;
    parent_module.add_submodule(thermodynamics)?;
    thermodynamics.add_class::<SWFJC>()?;
    Ok(())
}

/// The square-well freely-jointed chain (SWFJC) model thermodynamics.
#[pyclass]
#[derive(Copy, Clone)]
pub struct SWFJC
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

    /// The width of the well in units of nm.
    #[pyo3(get)]
    pub well_width: f64,

    /// The thermodynamic functions of the model in the isotensional ensemble.
    #[pyo3(get)]
    pub isotensional: super::isotensional::py::SWFJC
}

#[pymethods]
impl SWFJC
{
    #[new]
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, well_width: f64) -> Self
    {
        SWFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            well_width,
            isotensional: super::isotensional::py::SWFJC::init(number_of_links, link_length, hinge_mass, well_width)
        }
    }
}
