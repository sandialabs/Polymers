use pyo3::prelude::*;
use numpy::
{
    IntoPyArray,
    PyArrayDyn,
    PyReadonlyArrayDyn
};

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
    /// The expected force as a function of the applied end-to-end length and temperature,
    ///
    /// .. math::
    ///     f(\xi, T) = \frac{\partial \psi}{\partial\xi} = \frac{kT}{\xi} + \frac{kT}{\ell_b}\left(\frac{1}{2} - \frac{1}{N_b}\right)\frac{\sum_{s=0}^{s_\mathrm{max}}(-1)^s\binom{N_b}{s}\left(m - \frac{s}{N_b}\right)^{N_b - 3}}{\sum_{s=0}^{s_\mathrm{max}}(-1)^s\binom{N_b}{s}\left(m - \frac{s}{N_b}\right)^{N_b - 2}},
    ///
    /// where :math:`m\equiv(1 - \xi/N_b\ell_b)/2` and :math:`s_\mathrm{max}/N_b\leq m\leq (s_\mathrm{max}+1)/N_b`.
    ///
    /// Args:
    ///     end_to_end_length (numpy.ndarray): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The force :math:`f`.
    ///
    pub fn force<'py>(&self, py: Python<'py>, end_to_end_length: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).force(&end_to_end_length, &temperature)).into_pyarray(py)
    }
    /// The expected nondimensional force as a function of the applied nondimensional end-to-end length per link,
    ///
    /// .. math::
    ///     \eta(\gamma) = \frac{\partial\vartheta}{\partial\gamma} = \frac{1}{N_b\gamma} + \left(\frac{1}{2} - \frac{1}{N_b}\right)\frac{\sum_{s=0}^{s_\mathrm{max}}(-1)^s\binom{N_b}{s}\left(m - \frac{s}{N_b}\right)^{N_b - 3}}{\sum_{s=0}^{s_\mathrm{max}}(-1)^s\binom{N_b}{s}\left(m - \frac{s}{N_b}\right)^{N_b - 2}},
    ///
    /// where :math:`m\equiv(1 - \gamma)/2` and :math:`s_\mathrm{max}/N_b\leq m\leq (s_\mathrm{max}+1)/N_b`.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (numpy.ndarray): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///
    pub fn nondimensional_force<'py>(&self, py: Python<'py>, nondimensional_end_to_end_length_per_link: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_force(&nondimensional_end_to_end_length_per_link)).into_pyarray(py)
    }
    /// The Helmholtz free energy as a function of the applied end-to-end length and temperature,
    ///
    /// .. math::
    ///     \psi(\xi, T) = -kT\ln Q(\xi, T).
    ///
    /// Args:
    ///     end_to_end_length (numpy.ndarray): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The Helmholtz free energy :math:`\psi`.
    ///
    pub fn helmholtz_free_energy<'py>(&self, py: Python<'py>, end_to_end_length: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).helmholtz_free_energy(&end_to_end_length, &temperature)).into_pyarray(py)
    }
    /// The Helmholtz free energy per link as a function of the applied end-to-end length and temperature.
    ///
    /// Args:
    ///     end_to_end_length (numpy.ndarray): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The Helmholtz free energy per link :math:`\psi/N_b`.
    ///
    pub fn helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, end_to_end_length: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).helmholtz_free_energy_per_link(&end_to_end_length, &temperature)).into_pyarray(py)
    }
    /// The relative Helmholtz free energy as a function of the applied end-to-end length and temperature,
    ///
    /// .. math::
    ///     \Delta\psi(\xi, T) = kT\ln\left[\frac{P_\mathrm{eq}(0)}{P_\mathrm{eq}(\xi)}\right].
    ///
    /// Args:
    ///     end_to_end_length (numpy.ndarray): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The relative Helmholtz free energy :math:`\Delta\psi\equiv\psi(\xi,T)-\psi(0,T)`.
    ///
    pub fn relative_helmholtz_free_energy<'py>(&self, py: Python<'py>, end_to_end_length: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).relative_helmholtz_free_energy(&end_to_end_length, &temperature)).into_pyarray(py)
    }
    /// The relative Helmholtz free energy per link as a function of the applied end-to-end length and temperature.
    ///
    /// Args:
    ///     end_to_end_length (numpy.ndarray): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The relative Helmholtz free energy per link :math:`\Delta\psi/N_b`.
    ///
    pub fn relative_helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, end_to_end_length: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature)).into_pyarray(py)
    }
    /// The nondimensional Helmholtz free energy as a function of the applied nondimensional end-to-end length per link and temperature.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (numpy.ndarray): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional Helmholtz free energy :math:`\beta\psi=N_b\vartheta`.
    ///
    pub fn nondimensional_helmholtz_free_energy<'py>(&self, py: Python<'py>, nondimensional_end_to_end_length_per_link: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature)).into_pyarray(py)
    }
    /// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional end-to-end length per link and temperature.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (numpy.ndarray): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional Helmholtz free energy per link :math:`\vartheta\equiv\beta\psi/N_b`.
    ///
    pub fn nondimensional_helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_end_to_end_length_per_link: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature)).into_pyarray(py)
    }
    /// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional end-to-end length per link,
    ///
    /// .. math::
    ///     \beta\Delta\psi(\gamma) = \ln\left[\frac{\mathscr{P}_\mathrm{eq}(0)}{\mathscr{P}_\mathrm{eq}(\gamma)}\right].
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (numpy.ndarray): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Helmholtz free energy :math:`\beta\Delta\psi=N_b\Delta\vartheta`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy<'py>(&self, py: Python<'py>, nondimensional_end_to_end_length_per_link: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link)).into_pyarray(py)
    }
    /// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional end-to-end length per link,
    ///
    /// .. math::
    ///     \Delta\vartheta(\gamma) = \ln\left[\frac{\mathscr{P}_\mathrm{eq}(0)}{\mathscr{P}_\mathrm{eq}(\gamma)}\right]^{1/N_b}.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (numpy.ndarray): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Helmholtz free energy per link :math:`\Delta\vartheta\equiv\beta\Delta\psi/N_b`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_end_to_end_length_per_link: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link)).into_pyarray(py)
    }
    /// The nondimensional equilibrium probability density of nondimensional end-to-end vectors per link as a function of the nondimensional end-to-end length per link,
    ///
    /// .. math::
    ///     P_\mathrm{eq}(\xi) = \frac{e^{-\beta\psi(\xi, T)}}{4\pi\int e^{-\beta\psi(\xi', T)} \,{\xi'}{}^2 d\xi'} = \frac{1}{8\pi\ell_b^2\xi}\frac{N_b^{N_b - 2}}{(N_b - 2)!}\sum_{s=0}^{s_\mathrm{max}}(-1)^s\binom{N_b}{s}\left(m - \frac{s}{N_b}\right)^{N_b - 2},
    ///
    /// where :math:`m\equiv(1 - \xi/N_b\ell_b)/2` and :math:`s_\mathrm{max}/N_b\leq m\leq (s_\mathrm{max}+1)/N_b`.
    ///
    /// Args:
    ///     end_to_end_length (numpy.ndarray): The end-to-end length :math:`\xi`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The equilibrium probability density :math:`P_\mathrm{eq}`.
    ///
    pub fn equilibrium_distribution<'py>(&self, py: Python<'py>, end_to_end_length: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).equilibrium_distribution(&end_to_end_length)).into_pyarray(py)
    }
    /// The nondimensional equilibrium probability density of nondimensional end-to-end vectors per link as a function of the nondimensional end-to-end length per link,
    ///
    /// .. math::
    ///     \mathscr{P}_\mathrm{eq}(\gamma) = \frac{e^{-\Delta\vartheta(\gamma)}}{4\pi\int e^{-\Delta\vartheta(\gamma')} \,{\gamma'}{}^2 d\gamma'} = \frac{1}{8\pi\gamma}\frac{N_b^{N_b}}{(N_b - 2)!}\sum_{s=0}^{s_\mathrm{max}}(-1)^s\binom{N_b}{s}\left(m - \frac{s}{N_b}\right)^{N_b - 2},
    ///
    /// where :math:`m\equiv(1 - \gamma)/2` and :math:`s_\mathrm{max}/N_b\leq m\leq (s_\mathrm{max}+1)/N_b`.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (numpy.ndarray): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional equilibrium probability density :math:`\mathscr{P}_\mathrm{eq}\equiv (N_b\ell_b)^3 P_\mathrm{eq}`.
    ///
    pub fn nondimensional_equilibrium_distribution<'py>(&self, py: Python<'py>, nondimensional_end_to_end_length_per_link: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_equilibrium_distribution(&nondimensional_end_to_end_length_per_link)).into_pyarray(py)
    }
    /// The equilibrium probability density of end-to-end lengths as a function of the end-to-end length, given by :footcite:t:`treloar1949physics` as
    ///
    /// .. math::
    ///     g_\mathrm{eq}(\xi) = 4\pi\xi^2 P_\mathrm{eq}(\xi) = \frac{\xi}{2\ell_b^2}\frac{N_b^{N_b-2}}{(N_b - 2)!}\sum_{s=0}^{s_\mathrm{max}}(-1)^s\binom{N_b}{s}\left(m - \frac{s}{N_b}\right)^{N_b - 2},
    ///
    /// where :math:`m\equiv(1 - \xi/N_b\ell_b)/2` and :math:`s_\mathrm{max}/N_b\leq m\leq (s_\mathrm{max}+1)/N_b`.
    ///
    /// Args:
    ///     end_to_end_length (numpy.ndarray): The end-to-end length :math:`\xi`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The equilibrium probability density :math:`g_\mathrm{eq}`.
    ///
    pub fn equilibrium_radial_distribution<'py>(&self, py: Python<'py>, end_to_end_length: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).equilibrium_radial_distribution(&end_to_end_length)).into_pyarray(py)
    }
    /// The nondimensional equilibrium probability density of nondimensional end-to-end lengths per link as a function of the nondimensional end-to-end length per link,
    ///
    /// .. math::
    ///     \mathscr{g}_\mathrm{eq}(\gamma) = 4\pi\gamma^2 \mathscr{P}_\mathrm{eq}(\gamma) = \frac{\gamma}{2}\frac{N_b^{N_b}}{(N_b - 2)!}\sum_{s=0}^{s_\mathrm{max}}(-1)^s\binom{N_b}{s}\left(m - \frac{s}{N_b}\right)^{N_b - 2},
    ///
    /// where :math:`m\equiv(1 - \gamma)/2` and :math:`s_\mathrm{max}/N_b\leq m\leq (s_\mathrm{max}+1)/N_b`.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (numpy.ndarray): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional equilibrium probability density :math:`\mathscr{g}_\mathrm{eq}\equiv N_b\ell_b g_\mathrm{eq}`.
    ///
    pub fn nondimensional_equilibrium_radial_distribution<'py>(&self, py: Python<'py>, nondimensional_end_to_end_length_per_link: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_equilibrium_radial_distribution(&nondimensional_end_to_end_length_per_link)).into_pyarray(py)
    }
}
