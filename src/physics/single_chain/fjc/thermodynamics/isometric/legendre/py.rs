use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let legendre = PyModule::new(py, "legendre")?;
    parent_module.add_submodule(legendre)?;
    legendre.add_class::<FJC>()?;
    Ok(())
}

/// The freely-jointed chain (FJC) model thermodynamics in the isometric ensemble approximated using a Legendre transformation.
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
    /// The expected force as a function of the applied end-to-end length and temperature,
    ///
    /// .. math::
    ///     f(\xi, T) \sim \frac{kT}{\ell_b}\,\mathcal{L}^{-1}\left(\frac{\xi}{N_b\ell_b}\right) \quad \text{for } N_b\gg 1,
    ///
    /// where :math:`\mathcal{L}(x)=\coth(x)-1/x` is the Langevin function.
    ///
    /// Args:
    ///     end_to_end_length (float): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The force :math:`f`.
    ///
    pub fn force(&self, end_to_end_length_per_link: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).force(&end_to_end_length_per_link, &temperature))
    }
    /// The expected nondimensional force as a function of the applied nondimensional end-to-end length per link,
    ///
    /// .. math::
    ///     \eta(\gamma) \sim \mathcal{L}^{-1}(\gamma) \quad \text{for } N_b\gg 1,
    ///
    /// where :math:`\mathcal{L}(x)=\coth(x)-1/x` is the Langevin function.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (float): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     float: The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///
    pub fn nondimensional_force(&self, nondimensional_end_to_end_length_per_link: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_force(&nondimensional_end_to_end_length_per_link))
    }
    /// The helmholtz free energy as a function of the applied end-to-end length and temperature,
    ///
    /// .. math::
    ///     \psi(\xi, T) \sim \varphi\left[f(\xi, T)\right] + \xi f(\xi, T) \quad \text{for } N_b\gg 1,
    ///
    /// where :math:`f(\xi, T)` is given by the Legendre transformation approximation above.
    ///
    /// Args:
    ///     end_to_end_length (float): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The helmholtz free energy :math:`\psi`.
    ///
    pub fn helmholtz_free_energy(&self, end_to_end_length: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).helmholtz_free_energy(&end_to_end_length, &temperature))
    }
    /// The helmholtz free energy per link as a function of the applied end-to-end length and temperature.
    ///
    /// Args:
    ///     end_to_end_length (float): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The helmholtz free energy per link :math:`\psi/N_b`.
    ///
    pub fn helmholtz_free_energy_per_link(&self, end_to_end_length: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).helmholtz_free_energy_per_link(&end_to_end_length, &temperature))
    }
    /// The relative helmholtz free energy as a function of the applied end-to-end length and temperature.
    ///
    /// Args:
    ///     end_to_end_length (float): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The relative helmholtz free energy :math:`\Delta\psi\equiv\psi(\xi,T)-\psi(0,T)`.
    ///
    pub fn relative_helmholtz_free_energy(&self, end_to_end_length: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).relative_helmholtz_free_energy(&end_to_end_length, &temperature))
    }
    /// The relative helmholtz free energy per link as a function of the applied end-to-end length and temperature.
    ///
    /// Args:
    ///     end_to_end_length (float): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The relative helmholtz free energy per link :math:`\Delta\psi/N_b`.
    ///
    pub fn relative_helmholtz_free_energy_per_link(&self, end_to_end_length: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature))
    }
    /// The nondimensional helmholtz free energy as a function of the applied nondimensional end-to-end length per link and temperature.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (float): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The nondimensional helmholtz free energy :math:`\beta\psi=N_b\vartheta`.
    ///
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature))
    }
    /// The nondimensional helmholtz free energy per link as a function of the applied nondimensional end-to-end length per link and temperature,
    ///
    /// .. math::
    ///     \vartheta(\gamma, T) \sim \varphi\left[\mathcal{L}^{-1}(\gamma), T\right] + \gamma\mathcal{L}^{-1}(\gamma) \quad \text{for } N_b\gg 1,
    ///
    /// where :math:`\mathcal{L}(x)=\coth(x)-1/x` is the Langevin function.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (float): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The nondimensional helmholtz free energy per link :math:`\vartheta\equiv\beta\psi/N_b`.
    ///
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature))
    }
    /// The nondimensional relative helmholtz free energy as a function of the applied nondimensional end-to-end length per link.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (float): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     float: The nondimensional relative helmholtz free energy :math:`\beta\Delta\psi=N_b\Delta\vartheta`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link))
    }
    /// The nondimensional relative helmholtz free energy per link as a function of the applied nondimensional end-to-end length per link,
    ///
    /// .. math::
    ///     \Delta\vartheta(\gamma) \sim \gamma\mathcal{L}^{-1}(\gamma) + \ln\left\{\frac{\mathcal{L}^{-1}(\gamma)}{\sinh[\mathcal{L}^{-1}(\gamma)]}\right\} \quad \text{for } N_b\gg 1,
    ///
    /// where :math:`\mathcal{L}(x)=\coth(x)-1/x` is the Langevin function.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (float): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     float: The nondimensional relative helmholtz free energy per link :math:`\Delta\vartheta\equiv\beta\Delta\psi/N_b`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link))
    }
    /// The nondimensional equilibrium probability density of nondimensional end-to-end vectors per link as a function of the nondimensional end-to-end length per link,
    ///
    /// .. math::
    ///     P_\mathrm{eq}(\xi) = \frac{e^{-\beta\psi(\xi, T)}}{4\pi\int e^{-\beta\psi(\xi', T)} \,{\xi'}{}^2 d\xi'}.
    ///
    /// Args:
    ///     end_to_end_length (float): The end-to-end length :math:`\xi`.
    /// 
    /// Returns:
    ///     float: The equilibrium probability density :math:`P_\mathrm{eq}`.
    ///
    pub fn equilibrium_distribution(&self, end_to_end_length: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).equilibrium_distribution(&end_to_end_length))
    }
    /// The nondimensional equilibrium probability density of nondimensional end-to-end vectors per link as a function of the nondimensional end-to-end length per link,
    ///
    /// .. math::
    ///     \mathscr{P}_\mathrm{eq}(\gamma) = \frac{e^{-\Delta\vartheta(\gamma)}}{4\pi\int e^{-\Delta\vartheta(\gamma')} \,{\gamma'}{}^2 d\gamma'}.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (float): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     float: The nondimensional equilibrium probability density :math:`\mathscr{P}_\mathrm{eq}\equiv (N_b\ell_b)^3 P_\mathrm{eq}`.
    ///
    pub fn nondimensional_equilibrium_distribution(&self, nondimensional_end_to_end_length_per_link: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_equilibrium_distribution(&nondimensional_end_to_end_length_per_link))
    }
    /// The equilibrium probability density of end-to-end lengths as a function of the end-to-end length,
    ///
    /// .. math::
    ///     g_\mathrm{eq}(\xi) = 4\pi\xi^2 P_\mathrm{eq}(\xi).
    ///
    /// Args:
    ///     end_to_end_length (float): The end-to-end length :math:`\xi`.
    /// 
    /// Returns:
    ///     float: The equilibrium probability density :math:`g_\mathrm{eq}`.
    ///
    pub fn equilibrium_radial_distribution(&self, end_to_end_length: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).equilibrium_radial_distribution(&end_to_end_length))
    }
    /// The nondimensional equilibrium probability density of nondimensional end-to-end lengths per link as a function of the nondimensional end-to-end length per link,
    ///
    /// .. math::
    ///     \mathscr{g}_\mathrm{eq}(\gamma) = 4\pi\gamma^2 \mathscr{P}_\mathrm{eq}(\gamma).
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (float): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     float: The nondimensional equilibrium probability density :math:`\mathscr{g}_\mathrm{eq}\equiv N_b\ell_b g_\mathrm{eq}`.
    ///
    pub fn nondimensional_equilibrium_radial_distribution(&self, nondimensional_end_to_end_length_per_link: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_equilibrium_radial_distribution(&nondimensional_end_to_end_length_per_link))
    }
}
