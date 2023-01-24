use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let isotensional = PyModule::new(py, "isotensional")?;
    parent_module.add_submodule(isotensional)?;
    isotensional.add_class::<Ideal>()?;
    Ok(())
}

/// The ideal chain model thermodynamics in the isotensional ensemble.
#[pyclass]
#[derive(Copy, Clone)]
pub struct Ideal
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
impl Ideal
{
    #[new]
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        Ideal
        {
            hinge_mass,
            link_length,
            number_of_links
        }
    }
    /// The expected end-to-end length as a function of the applied force and temperature,
    ///
    /// .. math::
    ///     \xi(f, T) = -\frac{\partial\varphi}{\partial f} = \frac{N_b\ell_b^2f}{3kT}.
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The end-to-end length :math:`\xi`.
    ///
    pub fn end_to_end_length(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::Ideal::init(self.number_of_links, self.link_length, self.hinge_mass).end_to_end_length(&force, &temperature))
    }
    /// The expected end-to-end length per link as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The end-to-end length per link :math:`\xi/N_b=\ell_b\gamma`.
    ///
    pub fn end_to_end_length_per_link(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::Ideal::init(self.number_of_links, self.link_length, self.hinge_mass).end_to_end_length_per_link(&force, &temperature))
    }
    /// The expected nondimensional end-to-end length as a function of the applied nondimensional force.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional end-to-end length :math:`N_b\gamma=\xi/\ell_b`.
    ///
    pub fn nondimensional_end_to_end_length(&self, nondimensional_force: f64) -> PyResult<f64>
    {
        Ok(super::Ideal::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_end_to_end_length(&nondimensional_force))
    }
    /// The expected nondimensional end-to-end length per link as a function of the applied nondimensional force,
    ///
    /// .. math::
    ///     \gamma(\eta) = -\frac{\partial\varrho}{\partial\eta} = \frac{\eta}{3}.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    ///
    pub fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: f64) -> PyResult<f64>
    {
        Ok(super::Ideal::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_end_to_end_length_per_link(&nondimensional_force))
    }
    /// The gibbs free energy as a function of the applied force and temperature,
    ///
    /// .. math::
    ///     \varphi(f, T) = -kT\ln Z(f, T).
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The gibbs free energy :math:`\varphi`.
    ///
    pub fn gibbs_free_energy(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::Ideal::init(self.number_of_links, self.link_length, self.hinge_mass).gibbs_free_energy(&force, &temperature))
    }
    /// The gibbs free energy per link as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The gibbs free energy per link :math:`\varphi/N_b`.
    ///
    pub fn gibbs_free_energy_per_link(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::Ideal::init(self.number_of_links, self.link_length, self.hinge_mass).gibbs_free_energy_per_link(&force, &temperature))
    }
    /// The relative gibbs free energy as a function of the applied force and temperature,
    ///
    /// .. math::
    ///     \Delta\varphi(f, T) = -\frac{N_b\ell_b^2f^2}{6kT}.
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The relative gibbs free energy :math:`\Delta\varphi\equiv\varphi(f,T)-\varphi(0,T)`.
    ///
    pub fn relative_gibbs_free_energy(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::Ideal::init(self.number_of_links, self.link_length, self.hinge_mass).relative_gibbs_free_energy(&force, &temperature))
    }
    /// The relative gibbs free energy per link as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The relative gibbs free energy per link :math:`\Delta\varphi/N_b`.
    ///
    pub fn relative_gibbs_free_energy_per_link(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::Ideal::init(self.number_of_links, self.link_length, self.hinge_mass).relative_gibbs_free_energy_per_link(&force, &temperature))
    }
    /// The nondimensional gibbs free energy as a function of the applied nondimensional force and temperature.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional gibbs free energy :math:`\beta\varphi=N_b\varrho`.
    ///
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::Ideal::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_gibbs_free_energy(&nondimensional_force, &temperature))
    }
    /// The nondimensional gibbs free energy per link as a function of the applied nondimensional force and temperature.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional gibbs free energy per link :math:`\varrho\equiv\beta\varphi/N_b`.
    ///
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::Ideal::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature))
    }
    /// The nondimensional relative gibbs free energy as a function of the applied nondimensional force.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative gibbs free energy :math:`\beta\Delta\varphi=N_b\Delta\varrho`.
    ///
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: f64) -> PyResult<f64>
    {
        Ok(super::Ideal::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_relative_gibbs_free_energy(&nondimensional_force))
    }
    /// The nondimensional relative gibbs free energy per link as a function of the applied nondimensional force,
    ///
    /// .. math::
    ///     \Delta\varrho(\eta) = -\frac{\eta^2}{6}.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative gibbs free energy per link :math:`\Delta\varrho\equiv\beta\Delta\varphi/N_b`.
    ///
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: f64) -> PyResult<f64>
    {
        Ok(super::Ideal::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force))
    }
}
