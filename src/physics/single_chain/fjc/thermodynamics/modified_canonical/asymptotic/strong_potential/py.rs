use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let strong_potential = PyModule::new(py, "strong_potential")?;
    parent_module.add_submodule(strong_potential)?;
    strong_potential.add_class::<FJC>()?;
    Ok(())
}

/// The freely-jointed chain (FJC) model thermodynamics in the modified canonical ensemble approximated using an asymptotic approach valid for strong potentials.
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
    pub number_of_links: u8
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
            number_of_links
        }
    }
    /// The expected force as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (float): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The force :math:`f`.
    ///
    pub fn force(&self, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).force(&potential_distance, &potential_stiffness, &temperature))
    }
    /// The expected nondimensional force as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    ///
    /// Args:
    ///     nondimensional_potential_distance (float): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    /// 
    /// Returns:
    ///     float: The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///
    pub fn nondimensional_force(&self, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness))
    }
    /// The helmholtz free energy as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (float): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The helmholtz free energy :math:`\psi`.
    ///
    pub fn helmholtz_free_energy(&self, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature))
    }
    /// The helmholtz free energy per link as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (float): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The helmholtz free energy per link :math:`\psi/N_b`.
    ///
    pub fn helmholtz_free_energy_per_link(&self, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature))
    }
    /// The relative helmholtz free energy as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (float): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The relative helmholtz free energy :math:`\Delta\psi`.
    ///
    pub fn relative_helmholtz_free_energy(&self, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).relative_helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature))
    }
    /// The relative helmholtz free energy per link as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (float): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The relative helmholtz free energy per link :math:`\Delta\psi/N_b`.
    ///
    pub fn relative_helmholtz_free_energy_per_link(&self, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).relative_helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature))
    }
    /// The nondimensional helmholtz free energy as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    ///
    /// Args:
    ///     nondimensional_potential_distance (float): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The nondimensional helmholtz free energy :math:`\beta\psi=N_b\vartheta`.
    ///
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature))
    }
    /// The nondimensional helmholtz free energy per link as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    ///
    /// Args:
    ///     nondimensional_potential_distance (float): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The nondimensional helmholtz free energy per link :math:`\vartheta\equiv\beta\psi/N_b`.
    ///
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature))
    }
    /// The nondimensional relative helmholtz free energy as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    ///
    /// Args:
    ///     nondimensional_potential_distance (float): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    /// 
    /// Returns:
    ///     float: The nondimensional relative helmholtz free energy :math:`\beta\Delta\psi=N_b\Delta\vartheta`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_relative_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness))
    }
    /// The nondimensional relative helmholtz free energy per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    ///
    /// Args:
    ///     nondimensional_potential_distance (float): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    /// 
    /// Returns:
    ///     float: The nondimensional relative helmholtz free energy per link :math:`\Delta\vartheta\equiv\beta\Delta\psi/N_b`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness))
    }
}