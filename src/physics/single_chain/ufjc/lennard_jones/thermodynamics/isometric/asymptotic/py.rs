use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let asymptotic = PyModule::new(py, "asymptotic")?;
    super::reduced::py::register_module(py, asymptotic)?;
    super::legendre::py::register_module(py, asymptotic)?;
    parent_module.add_submodule(asymptotic)?;
    asymptotic.add_class::<LENNARDJONESFJC>()?;
    Ok(())
}

/// The Lennard-Jones link potential freely-jointed chain (Lennard-Jones-FJC) model thermodynamics in the isometric ensemble approximated using an asymptotic approach.
#[pyclass]
#[derive(Copy, Clone)]
pub struct LENNARDJONESFJC
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

    /// The thermodynamic functions of the model in the isometric ensemble approximated using a reduced asymptotic approach.
    #[pyo3(get)]
    pub reduced: super::reduced::py::LENNARDJONESFJC,

    /// The thermodynamic functions of the model in the isometric ensemble approximated using an asymptotic approach and a Legendre transformation.
    #[pyo3(get)]
    pub legendre: super::legendre::py::LENNARDJONESFJC
}

#[pymethods]
impl LENNARDJONESFJC
{
    #[new]
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        LENNARDJONESFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            reduced: super::reduced::py::LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness),
            legendre: super::legendre::py::LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness)
        }
    }
}
