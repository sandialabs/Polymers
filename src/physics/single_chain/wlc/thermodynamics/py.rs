use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let thermodynamics = PyModule::new(py, "thermodynamics")?;
    super::isometric::py::register_module(py, thermodynamics)?;
    parent_module.add_submodule(thermodynamics)?;
    thermodynamics.add_class::<WLC>()?;
    Ok(())
}

/// The worm-like chain (WLC) model thermodynamics.
#[pyclass]
#[derive(Copy, Clone)]
pub struct WLC
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

    /// The persistance length of the chain in units of nm.
    #[pyo3(get)]
    pub persistance_length: f64,

    /// The thermodynamic functions of the model in the isometric ensemble.
    #[pyo3(get)]
    pub isometric: super::isometric::py::WLC
}

#[pymethods]
impl WLC
{
    #[new]
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, persistance_length: f64) -> Self
    {
        WLC
        {
            hinge_mass,
            link_length,
            number_of_links,
            persistance_length,
            isometric: super::isometric::py::WLC::init(number_of_links, link_length, hinge_mass, persistance_length)
        }
    }
}
