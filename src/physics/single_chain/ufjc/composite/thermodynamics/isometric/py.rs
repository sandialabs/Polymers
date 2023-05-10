use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let isometric = PyModule::new(py, "isometric")?;
    super::legendre::py::register_module(py, isometric)?;
    parent_module.add_submodule(isometric)?;
    isometric.add_class::<CUFJC>()?;
    Ok(())
}

/// The composite uFJC (CuFJC) model thermodynamics in the isometric ensemble.
#[pyclass]
#[derive(Copy, Clone)]
pub struct CUFJC
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

    /// The number of bonds in each link.
    #[pyo3(get)]
    pub number_of_bonds: u8,

    /// The stiffness of each bond in units of J/(mol⋅nm^2).
    #[pyo3(get)]
    pub bond_stiffness: f64,

    /// The energy of each bond in units of J/mol.
    #[pyo3(get)]
    pub bond_energy: f64,

    /// The scission energy of each bond in units of J/mol.
    #[pyo3(get)]
    pub bond_scission_energy: f64,

    /// The attempt frequency of each bond in units of 1/ns.
    #[pyo3(get)]
    pub bond_attempt_frequency: f64,
    
    /// The thermodynamic functions of the model in the isometric ensemble approximated using a Legendre transformation.
    #[pyo3(get)]
    pub legendre: super::legendre::py::CUFJC
}

#[pymethods]
impl CUFJC
{
    #[new]
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, number_of_bonds: u8, bond_stiffness: f64, bond_energy: f64, bond_scission_energy: f64, bond_attempt_frequency: f64) -> Self
    {
        CUFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_bonds,
            bond_stiffness,
            bond_energy,
            bond_scission_energy,
            bond_attempt_frequency,
            legendre: super::legendre::py::CUFJC::init(number_of_links, link_length, hinge_mass, number_of_bonds, bond_stiffness, bond_energy, bond_scission_energy, bond_attempt_frequency)
        }
    }
}
