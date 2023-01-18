use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let isometric = PyModule::new(py, "isometric")?;
    super::legendre::py::register_module(py, isometric)?;
    parent_module.add_submodule(isometric)?;
    isometric.add_class::<FJC>()?;
    Ok(())
}

/// The freely-jointed chain (FJC) model thermodynamics in the isometric ensemble.
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
    pub number_of_links: u8,

    /// The thermodynamic functions of the model in the isometric ensemble approximated using a Legendre transformation.
    /// 
    /// Example: 
    ///     Plot the nondimensional force :math:`\eta` as a function of the nondimensional end-to-end length per link :math:`\gamma` for an increasing number of links :math:`N_b` and compare with the approximation obtained when using a Legendre transformation:
    /// 
    ///     .. plot::
    /// 
    ///         >>> import numpy as np
    ///         >>> import matplotlib.pyplot as plt
    ///         >>> from polymers import physics
    ///         >>> FJC = physics.single_chain.fjc.thermodynamics.isometric.FJC
    ///         >>> gamma = np.linspace(0, 1, 100)[1:-1]
    ///         >>> eta = np.zeros(len(gamma))
    ///         >>> eta_legendre = np.zeros(len(gamma))
    ///         >>> for N_b in [4, 8, 16, 32]:
    ///         ...     fjc = FJC(N_b, 1, 1)
    ///         ...     for i, gamma_i in enumerate(gamma):
    ///         ...         eta[i] = fjc.nondimensional_force(gamma_i)
    ///         ...     plt.plot(gamma, eta, label=r'$N_b=$' + str(N_b))
    ///         >>> for i, gamma_i in enumerate(gamma):
    ///         ...     eta_legendre[i] = fjc.legendre.nondimensional_force(gamma_i)
    ///         >>> plt.plot(gamma, eta_legendre, 'k--', label='legendre')
    ///         >>> plt.legend()
    ///         >>> plt.xlim([0, 1])
    ///         >>> plt.ylim([0, 10])
    ///         >>> plt.xlabel(r'$\gamma$')
    ///         >>> plt.ylabel(r'$\eta$')
    ///         >>> plt.show()
    ///
    #[pyo3(get)]
    pub legendre: super::legendre::py::FJC
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
            number_of_links,
            legendre: super::legendre::py::FJC::init(number_of_links, link_length, hinge_mass)
        }
    }
    /// The force as a function of the end-to-end length and temperature,
    ///
    /// .. math::
    ///     f(\xi, T) = \frac{\partial \psi}{\partial\xi} = ?.
    ///
    /// Args:
    ///     end_to_end_length (float): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     float: The force :math:`f`.
    ///
    pub fn force(&self, end_to_end_length: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).force(&end_to_end_length, &temperature))
    }
    /// The nondimensional force as a function of the nondimensional end-to-end length per link,
    ///
    /// .. math::
    ///     \eta(\gamma) = \frac{\partial\vartheta}{\partial\gamma} = ?.
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
    /// The equilibrium probability density of end-to-end lengths as a function of the end-to-end length, given by :cite:t:`treloar1949physics` as
    ///
    /// .. math::
    ///     g_\mathrm{eq}(\xi) = 4\pi\xi^2 P_\mathrm{eq}(\xi) = \frac{\xi}{2\ell_b^2}\frac{N_b^{N_b-2}}{(N_b - 2)!}\sum_{s=0}^{k}(-1)^s\binom{N_b}{s}\left(m - \frac{s}{N_b}\right)^{N_b - 2},
    ///
    /// where :math:`m\equiv(1 - \xi/N_b\ell_b)/2` and :math:`k/N_b\leq m\leq (k+1)/N_b`.
    ///
    /// Args:
    ///     end_to_end_length (float): The end-to-end length :math:`\xi`.
    /// 
    /// Returns:
    ///     float: The equilibrium probability density :math:`g_\mathrm{eq}`.
    ///
    pub fn equilibrium_radial_distribution(&self, end_to_end_length: f64) -> PyResult<f64>
    {
        Ok(super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_equilibrium_radial_distribution(&end_to_end_length))
    }
}
