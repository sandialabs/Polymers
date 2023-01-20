use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let legendre = PyModule::new(py, "legendre")?;
    parent_module.add_submodule(legendre)?;
    legendre.add_class::<SWFJC>()?;
    Ok(())
}

/// The square-well freely-jointed chain (SWFJC) model thermodynamics in the isotensional ensemble approximated using a Legendre transformation.
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
    pub well_width: f64
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
            well_width
        }
    }
    /// The helmholtz free energy as a function of the applied force and temperature,
    ///
    /// .. math::
    ///     \psi(f, T) \sim \varphi(f, T) + f \xi(f, T) \quad \text{for } N_b\gg 1.
    ///
    /// Args:
    ///     force (float): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The helmholtz free energy :math:`\psi`.
    ///
    pub fn helmholtz_free_energy(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::SWFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.well_width).helmholtz_free_energy(&force, &temperature))
    }
    /// The helmholtz free energy per link as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (float): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The helmholtz free energy per link :math:`\psi/N_b`.
    ///
    pub fn helmholtz_free_energy_per_link(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::SWFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.well_width).helmholtz_free_energy_per_link(&force, &temperature))
    }
    /// The relative helmholtz free energy as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (float): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The relative helmholtz free energy :math:`\Delta\psi\equiv\psi(f,T)-\psi(0,T)`.
    ///
    pub fn relative_helmholtz_free_energy(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::SWFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.well_width).relative_helmholtz_free_energy(&force, &temperature))
    }
    /// The relative helmholtz free energy per link as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (float): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The relative helmholtz free energy per link :math:`\Delta\psi/N_b`.
    ///
    pub fn relative_helmholtz_free_energy_per_link(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::SWFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.well_width).relative_helmholtz_free_energy_per_link(&force, &temperature))
    }
    /// The nondimensional helmholtz free energy as a function of the applied nondimensional force and temperature.
    ///
    /// Args:
    ///     nondimensional_force (float): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The nondimensional helmholtz free energy :math:`\beta\psi=N_b\vartheta`.
    ///
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::SWFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.well_width).nondimensional_helmholtz_free_energy(&nondimensional_force, &temperature))
    }
    /// The nondimensional helmholtz free energy per link as a function of the applied nondimensional force and temperature.
    ///
    /// Args:
    ///     nondimensional_force (float): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The nondimensional helmholtz free energy per link :math:`\vartheta\equiv\beta\psi/N_b`.
    ///
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::SWFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.well_width).nondimensional_helmholtz_free_energy_per_link(&nondimensional_force, &temperature))
    }
    /// The nondimensional relative helmholtz free energy as a function of the applied nondimensional force.
    ///
    /// Args:
    ///     nondimensional_force (float): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    /// 
    /// Returns:
    ///     float: The nondimensional relative helmholtz free energy :math:`\beta\Delta\psi=N_b\Delta\vartheta`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_force: f64) -> PyResult<f64>
    {
        Ok(super::SWFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.well_width).nondimensional_relative_helmholtz_free_energy(&nondimensional_force))
    }
    /// The nondimensional relative helmholtz free energy per link as a function of the applied nondimensional force.
    ///
    /// Args:
    ///     nondimensional_force (float): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    /// 
    /// Returns:
    ///     float: The nondimensional relative helmholtz free energy per link :math:`\Delta\vartheta\equiv\beta\Delta\psi/N_b`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_force: f64) -> PyResult<f64>
    {
        Ok(super::SWFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.well_width).nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force))
    }
}
