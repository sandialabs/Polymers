use pyo3::prelude::*;
use numpy::
{
    IntoPyArray,
    PyArrayDyn,
    PyReadonlyArrayDyn
};

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let isotensional = PyModule::new(py, "isotensional")?;
    super::legendre::py::register_module(py, isotensional)?;
    parent_module.add_submodule(isotensional)?;
    isotensional.add_class::<SWFJC>()?;
    Ok(())
}

/// The square-well freely-jointed chain (SWFJC) model thermodynamics in the isotensional ensemble.
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
    pub well_width: f64,

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using a Legendre transformation.
    #[pyo3(get)]
    pub legendre: super::legendre::py::SWFJC
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
            well_width,
            legendre: super::legendre::py::SWFJC::init(number_of_links, link_length, hinge_mass, well_width)
        }
    }
    /// The expected end-to-end length as a function of the applied force and temperature,
    ///
    /// .. math::
    ///     \xi(f, T) = -\frac{\partial\varphi}{\partial f}.
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The end-to-end length :math:`\xi`.
    ///
    pub fn end_to_end_length<'py>(&self, py: Python<'py>, force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        force.as_array().mapv(|force: f64| super::end_to_end_length(&self.number_of_links, &self.link_length, &self.well_width, &force, &temperature)).into_pyarray(py)
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
    pub fn end_to_end_length_per_link<'py>(&self, py: Python<'py>, force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        force.as_array().mapv(|force: f64| super::end_to_end_length_per_link(&self.link_length, &self.well_width, &force, &temperature)).into_pyarray(py)
    }
    /// The expected nondimensional end-to-end length as a function of the applied nondimensional force.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional end-to-end length :math:`N_b\gamma=\xi/\ell_b`.
    ///
    pub fn nondimensional_end_to_end_length<'py>(&self, py: Python<'py>, nondimensional_force: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_end_to_end_length(&self.number_of_links, &self.link_length, &self.well_width, &nondimensional_force)).into_pyarray(py)
    }
    /// The expected nondimensional end-to-end length per link as a function of the applied nondimensional force, given by :footcite:t:`buche2022freely` as
    ///
    /// .. math::
    ///     \gamma(\eta) = -\frac{\partial\varrho}{\partial\eta} = \frac{\varsigma^2\eta\sinh(\varsigma\eta) - \eta\sinh(\eta)}{y(\eta,\varsigma) - y(\eta,1)},
    ///
    /// where :math:`\varsigma\equiv 1+w/\ell_b` is the nondimensional well width parameter,
    /// and where :math:`y(\eta,\varsigma)\equiv\varsigma\eta\cosh(\varsigma\eta)-\sinh(\varsigma\eta)`.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    ///
    pub fn nondimensional_end_to_end_length_per_link<'py>(&self, py: Python<'py>, nondimensional_force: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_end_to_end_length_per_link(&self.link_length, &self.well_width, &nondimensional_force)).into_pyarray(py)
    }
    /// The Gibbs free energy as a function of the applied force and temperature,
    ///
    /// .. math::
    ///     \varphi(f, T) = -kT\ln Z(f, T).
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The Gibbs free energy :math:`\varphi`.
    ///
    pub fn gibbs_free_energy<'py>(&self, py: Python<'py>, force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        force.as_array().mapv(|force: f64| super::gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.well_width, &force, &temperature)).into_pyarray(py)
    }
    /// The Gibbs free energy per link as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The Gibbs free energy per link :math:`\varphi/N_b`.
    ///
    pub fn gibbs_free_energy_per_link<'py>(&self, py: Python<'py>, force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        force.as_array().mapv(|force: f64| super::gibbs_free_energy_per_link(&self.link_length, &self.hinge_mass, &self.well_width, &force, &temperature)).into_pyarray(py)
    }
    /// The relative Gibbs free energy as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The relative Gibbs free energy :math:`\Delta\varphi\equiv\varphi(f,T)-\varphi(0,T)`.
    ///
    pub fn relative_gibbs_free_energy<'py>(&self, py: Python<'py>, force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        force.as_array().mapv(|force: f64| super::relative_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.well_width, &force, &temperature)).into_pyarray(py)
    }
    /// The relative Gibbs free energy per link as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The relative Gibbs free energy per link :math:`\Delta\varphi/N_b`.
    ///
    pub fn relative_gibbs_free_energy_per_link<'py>(&self, py: Python<'py>, force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        force.as_array().mapv(|force: f64| super::relative_gibbs_free_energy_per_link(&self.link_length, &self.well_width, &force, &temperature)).into_pyarray(py)
    }
    /// The nondimensional Gibbs free energy as a function of the applied nondimensional force and temperature.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional Gibbs free energy :math:`\beta\varphi=N_b\varrho`.
    ///
    pub fn nondimensional_gibbs_free_energy<'py>(&self, py: Python<'py>, nondimensional_force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.well_width, &nondimensional_force, &temperature)).into_pyarray(py)
    }
    /// The nondimensional Gibbs free energy per link as a function of the applied nondimensional force and temperature.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional Gibbs free energy per link :math:`\varrho\equiv\beta\varphi/N_b`.
    ///
    pub fn nondimensional_gibbs_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_gibbs_free_energy_per_link(&self.link_length, &self.hinge_mass, &self.well_width, &nondimensional_force, &temperature)).into_pyarray(py)
    }
    /// The nondimensional relative Gibbs free energy as a function of the applied nondimensional force.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Gibbs free energy :math:`\beta\Delta\varphi=N_b\Delta\varrho`.
    ///
    pub fn nondimensional_relative_gibbs_free_energy<'py>(&self, py: Python<'py>, nondimensional_force: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_relative_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.well_width, &nondimensional_force)).into_pyarray(py)
    }
    /// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional force, given by :footcite:t:`buche2022freely` as
    ///
    /// .. math::
    ///     \Delta\varrho(\eta) = 3\ln(\eta) - \ln\left[y(\eta,\varsigma) - y(\eta,1)\right],
    ///
    /// where :math:`\varsigma\equiv 1+w/\ell_b` is the nondimensional well width parameter,
    /// and where :math:`y(\eta,\varsigma)\equiv\varsigma\eta\cosh(\varsigma\eta)-\sinh(\varsigma\eta)`.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Gibbs free energy per link :math:`\Delta\varrho\equiv\beta\Delta\varphi/N_b`.
    ///
    pub fn nondimensional_relative_gibbs_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_force: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_relative_gibbs_free_energy_per_link(&self.link_length, &self.well_width, &nondimensional_force)).into_pyarray(py)
    }
}
