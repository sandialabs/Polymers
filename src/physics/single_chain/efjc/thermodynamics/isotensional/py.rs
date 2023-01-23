use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let isotensional = PyModule::new(py, "isotensional")?;
    super::asymptotic::py::register_module(py, isotensional)?;
    super::legendre::py::register_module(py, isotensional)?;
    parent_module.add_submodule(isotensional)?;
    isotensional.add_class::<EFJC>()?;
    Ok(())
}

/// The extensible freely-jointed chain (EFJC) model thermodynamics in the isotensional ensemble.
#[pyclass]
#[derive(Copy, Clone)]
pub struct EFJC
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

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using an asymptotic approach.
    #[pyo3(get)]
    pub asymptotic: super::asymptotic::py::EFJC,

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using a Legendre transformation.
    #[pyo3(get)]
    pub legendre: super::legendre::py::EFJC
}

#[pymethods]
impl EFJC
{
    #[new]
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        EFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            asymptotic: super::asymptotic::py::EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness),
            legendre: super::legendre::py::EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness)
        }
    }
    /// The expected end-to-end length as a function of the applied force and temperature,
    ///
    /// .. math::
    ///     \xi(f, T) = -\frac{\partial\varphi}{\partial f}.
    ///
    /// Args:
    ///     force (float): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The end-to-end length :math:`\xi`.
    ///
    pub fn end_to_end_length(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::EFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.link_stiffness).end_to_end_length(&force, &temperature))
    }
    /// The expected end-to-end length per link as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (float): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The end-to-end length per link :math:`\xi/N_b=\ell_b\gamma`.
    ///
    pub fn end_to_end_length_per_link(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::EFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.link_stiffness).end_to_end_length_per_link(&force, &temperature))
    }
    /// The expected nondimensional end-to-end length as a function of the applied nondimensional force.
    ///
    /// Args:
    ///     nondimensional_force (float): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The nondimensional end-to-end length :math:`N_b\gamma=\xi/\ell_b`.
    ///
    pub fn nondimensional_end_to_end_length(&self, nondimensional_force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::EFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.link_stiffness).nondimensional_end_to_end_length(&nondimensional_force, &temperature))
    }
    /// The expected nondimensional end-to-end length per link as a function of the applied nondimensional force, calculated from :cite:t:`balabaev2009extension,buche2022freely` as
    ///
    /// .. math::
    ///     \gamma(\eta) = -\frac{\partial\varrho}{\partial\eta} = \mathcal{L}(\eta) + \frac{\eta}{\kappa}\left[1 + \frac{1 - \mathcal{L}(\eta)\coth(\eta)}{1 + (\eta/\kappa)\coth(\eta)}\right] + \frac{\partial}{\partial\eta}\,\ln\left[1+g(\eta)\right],
    ///
    /// where :math:`\mathcal{L}(x)=\coth(x)-1/x` is the Langevin function, and :math:`g(\eta)` is defined as
    ///
    /// .. math::
    ///     g(\eta) \equiv \frac{e^{\eta}\left(\frac{\eta}{\kappa} + 1\right) \,\mathrm{erf}\left(\frac{\eta+\kappa}{\sqrt{2\kappa}}\right) - e^{-\eta}\left(\frac{\eta}{\kappa} - 1\right) \,\mathrm{erf}\left(\frac{\eta-\kappa}{\sqrt{2\kappa}}\right)}{4\sinh(\eta)\left[1 + (\eta/\kappa)\coth(\eta)\right]} - \frac{1}{2} .
    ///
    /// Args:
    ///     nondimensional_force (float): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    ///
    pub fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::EFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.link_stiffness).nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature))
    }
    /// The gibbs free energy as a function of the applied force and temperature,
    ///
    /// .. math::
    ///     \varphi(f, T) = -kT\ln Z(f, T).
    ///
    /// Args:
    ///     force (float): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The gibbs free energy :math:`\varphi`.
    ///
    pub fn gibbs_free_energy(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::EFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.link_stiffness).gibbs_free_energy(&force, &temperature))
    }
    /// The gibbs free energy per link as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (float): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The gibbs free energy per link :math:`\varphi/N_b`.
    ///
    pub fn gibbs_free_energy_per_link(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::EFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.link_stiffness).gibbs_free_energy_per_link(&force, &temperature))
    }
    /// The relative gibbs free energy as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (float): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The relative gibbs free energy :math:`\Delta\varphi\equiv\varphi(f,T)-\varphi(0,T)`.
    ///
    pub fn relative_gibbs_free_energy(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::EFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.link_stiffness).relative_gibbs_free_energy(&force, &temperature))
    }
    /// The relative gibbs free energy per link as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (float): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The relative gibbs free energy per link :math:`\Delta\varphi/N_b`.
    ///
    pub fn relative_gibbs_free_energy_per_link(&self, force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::EFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.link_stiffness).relative_gibbs_free_energy_per_link(&force, &temperature))
    }
    /// The nondimensional gibbs free energy as a function of the applied nondimensional force and temperature.
    ///
    /// Args:
    ///     nondimensional_force (float): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The nondimensional gibbs free energy :math:`\beta\varphi=N_b\varrho`.
    ///
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::EFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.link_stiffness).nondimensional_gibbs_free_energy(&nondimensional_force, &temperature))
    }
    /// The nondimensional gibbs free energy per link as a function of the applied nondimensional force and temperature.
    ///
    /// Args:
    ///     nondimensional_force (float): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The nondimensional gibbs free energy per link :math:`\varrho\equiv\beta\varphi/N_b`.
    ///
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::EFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.link_stiffness).nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature))
    }
    /// The nondimensional relative gibbs free energy as a function of the applied nondimensional force.
    ///
    /// Args:
    ///     nondimensional_force (float): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The nondimensional relative gibbs free energy :math:`\beta\Delta\varphi=N_b\Delta\varrho`.
    ///
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::EFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.link_stiffness).nondimensional_relative_gibbs_free_energy(&nondimensional_force, &temperature))
    }
    /// The nondimensional relative gibbs free energy per link as a function of the applied nondimensional force, given by :cite:t:`buche2022freely` as
    ///
    /// .. math::
    ///     \Delta\varrho(\eta) = -\frac{\eta^2}{2\kappa} - \ln\left[w^+(\eta) + w^-(\eta)\right],
    ///
    /// where the functions :math:`w^+(\eta)` and :math:`w^-(\eta)` are defined as as
    ///
    /// .. math::
    ///     w^\pm(\eta) \equiv e^{\pm\eta}\left(\frac{1}{\kappa} \pm \frac{1}{\eta}\right) \left[1 \pm \mathrm{erf}\left(\frac{\eta\pm\kappa}{\sqrt{2\kappa}}\right).\right]
    ///
    /// Args:
    ///     nondimensional_force (float): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The nondimensional relative gibbs free energy per link :math:`\Delta\varrho\equiv\beta\Delta\varphi/N_b`.
    ///
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::EFJC::init(self.number_of_links, self.link_length, self.hinge_mass, self.link_stiffness).nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force, &temperature))
    }
}
