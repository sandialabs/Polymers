use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let weak_potential = PyModule::new(py, "weak_potential")?;
    parent_module.add_submodule(weak_potential)?;
    weak_potential.add_class::<FJC>()?;
    Ok(())
}

/// The freely-jointed chain (FJC) model thermodynamics in the modified canonical ensemble approximated using an asymptotic approach valid for weak potentials.
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
    /// The expected force as a function of the applied potential distance, and potential stiffness.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    /// 
    /// Returns:
    ///     numpy.ndarray: The force :math:`f`.
    ///
    pub fn force(&self, potential_distance: f64, potential_stiffness: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).force(&potential_distance, &potential_stiffness))
    }
    /// The expected nondimensional force as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///
    pub fn nondimensional_force(&self, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness))
    }
    /// The expected end-to-end length as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The end-to-end length :math:`\xi`.
    ///
    pub fn end_to_end_length(&self, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).end_to_end_length(&potential_distance, &potential_stiffness, &temperature))
    }
    /// The expected end-to-end length per link as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The end-to-end length per link :math:`\xi/N_b=\ell_b\gamma`.
    ///
    pub fn end_to_end_length_per_link(&self, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).end_to_end_length_per_link(&potential_distance, &potential_stiffness, &temperature))
    }
    /// The expected nondimensional end-to-end length as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional end-to-end length :math:`N_b\gamma=\xi/\ell_b`.
    ///
    pub fn nondimensional_end_to_end_length(&self, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_end_to_end_length(&nondimensional_potential_distance, &nondimensional_potential_stiffness))
    }
    /// The expected nondimensional end-to-end length per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional end-to-end length :math:`\gamma\equiv\xi/N_b\ell_b`.
    ///
    pub fn nondimensional_end_to_end_length_per_link(&self, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_end_to_end_length_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness))
    }
    /// The gibbs free energy as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The gibbs free energy :math:`\varphi`.
    ///
    pub fn gibbs_free_energy(&self, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature))
    }
    /// The gibbs free energy epr link as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The gibbs free energy per link :math:`\varphi/N_b`.
    ///
    pub fn gibbs_free_energy_per_link(&self, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature))
    }
    /// The relative gibbs free energy as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The relative gibbs free energy :math:`\Delta\varphi`.
    ///
    pub fn relative_gibbs_free_energy(&self, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).relative_gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature))
    }
    /// The relative gibbs free energy per link as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The relative gibbs free energy per link :math:`\Delta\varphi/N_b`.
    ///
    pub fn relative_gibbs_free_energy_per_link(&self, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).relative_gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature))
    }
    /// The nondimensional gibbs free energy as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional gibbs free energy :math:`\beta\varphi=N_b\varrho`.
    ///
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature))
    }
    /// The nondimensional gibbs free energy per link as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional gibbs free energy per link :math:`\varrho\equiv\beta\varphi/N_b`.
    ///
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature))
    }
    /// The nondimensional relative gibbs free energy as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative gibbs free energy :math:`\beta\Delta\varphi=N_b\Delta\varrho`.
    ///
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness))
    }
    /// The nondimensional relative gibbs free energy per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative gibbs free energy per link :math:`\Delta\varrho\equiv\beta\Delta\varphi/N_b`.
    ///
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness))
    }
}