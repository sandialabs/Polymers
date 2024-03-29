use pyo3::prelude::*;
use numpy::
{
    IntoPyArray,
    PyArrayDyn,
    PyReadonlyArrayDyn
};

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let modified_canonical = PyModule::new(py, "modified_canonical")?;
    super::asymptotic::py::register_module(py, modified_canonical)?;
    parent_module.add_submodule(modified_canonical)?;
    modified_canonical.add_class::<FJC>()?;
    Ok(())
}

/// The freely-jointed chain (FJC) model thermodynamics in the modified canonical ensemble.
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

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using an asymptotic approach.
    #[pyo3(get)]
    pub asymptotic: super::asymptotic::py::FJC
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
            asymptotic: super::asymptotic::py::FJC::init(number_of_links, link_length, hinge_mass)
        }
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
    pub fn end_to_end_length<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::end_to_end_length(&self.number_of_links, &self.link_length, &potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
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
    pub fn end_to_end_length_per_link<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::end_to_end_length_per_link(&self.number_of_links, &self.link_length, &potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
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
    pub fn nondimensional_end_to_end_length<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::nondimensional_end_to_end_length(&self.number_of_links, &nondimensional_potential_distance, &nondimensional_potential_stiffness)).into_pyarray(py)
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
    pub fn nondimensional_end_to_end_length_per_link<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::nondimensional_end_to_end_length_per_link(&self.number_of_links, &nondimensional_potential_distance, &nondimensional_potential_stiffness)).into_pyarray(py)
    }
    /// The expected force as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The force :math:`f`.
    ///
    pub fn force<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::force(&self.number_of_links, &self.link_length, &potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
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
    pub fn nondimensional_force<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::nondimensional_force(&self.number_of_links, &nondimensional_potential_distance, &nondimensional_potential_stiffness)).into_pyarray(py)
    }
    /// The Helmholtz free energy as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The Helmholtz free energy :math:`\psi`.
    ///
    pub fn helmholtz_free_energy<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The Helmholtz free energy per link as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The Helmholtz free energy per link :math:`\psi/N_b`.
    ///
    pub fn helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The relative Helmholtz free energy as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The relative Helmholtz free energy :math:`\Delta\psi`.
    ///
    pub fn relative_helmholtz_free_energy<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::relative_helmholtz_free_energy(&self.number_of_links, &self.link_length, &potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The relative Helmholtz free energy per link as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The relative Helmholtz free energy per link :math:`\Delta\psi/N_b`.
    ///
    pub fn relative_helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::relative_helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The nondimensional Helmholtz free energy as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The nondimensional Helmholtz free energy :math:`\beta\psi=N_b\vartheta`.
    ///
    pub fn nondimensional_helmholtz_free_energy<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::nondimensional_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The nondimensional Helmholtz free energy per link :math:`\vartheta\equiv\beta\psi/N_b`.
    ///
    pub fn nondimensional_helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::nondimensional_helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Helmholtz free energy :math:`\beta\Delta\psi=N_b\Delta\vartheta`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::nondimensional_relative_helmholtz_free_energy(&self.number_of_links, &nondimensional_potential_distance, &nondimensional_potential_stiffness)).into_pyarray(py)
    }
    /// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Helmholtz free energy per link :math:`\Delta\vartheta\equiv\beta\Delta\psi/N_b`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::nondimensional_relative_helmholtz_free_energy_per_link(&self.number_of_links, &nondimensional_potential_distance, &nondimensional_potential_stiffness)).into_pyarray(py)
    }
    /// The Gibbs free energy as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The Gibbs free energy :math:`\varphi`.
    ///
    pub fn gibbs_free_energy<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The Gibbs free energy epr link as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The Gibbs free energy per link :math:`\varphi/N_b`.
    ///
    pub fn gibbs_free_energy_per_link<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The relative Gibbs free energy as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The relative Gibbs free energy :math:`\Delta\varphi`.
    ///
    pub fn relative_gibbs_free_energy<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::relative_gibbs_free_energy(&self.number_of_links, &self.link_length, &potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The relative Gibbs free energy per link as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The relative Gibbs free energy per link :math:`\Delta\varphi/N_b`.
    ///
    pub fn relative_gibbs_free_energy_per_link<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::relative_gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The nondimensional Gibbs free energy as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The nondimensional Gibbs free energy :math:`\beta\varphi=N_b\varrho`.
    ///
    pub fn nondimensional_gibbs_free_energy<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::nondimensional_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The nondimensional Gibbs free energy per link as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The nondimensional Gibbs free energy per link :math:`\varrho\equiv\beta\varphi/N_b`.
    ///
    pub fn nondimensional_gibbs_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::nondimensional_gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The nondimensional relative Gibbs free energy as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Gibbs free energy :math:`\beta\Delta\varphi=N_b\Delta\varrho`.
    ///
    pub fn nondimensional_relative_gibbs_free_energy<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::nondimensional_relative_gibbs_free_energy(&self.number_of_links, &nondimensional_potential_distance, &nondimensional_potential_stiffness)).into_pyarray(py)
    }
    /// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Gibbs free energy per link :math:`\Delta\varrho\equiv\beta\Delta\varphi/N_b`.
    ///
    pub fn nondimensional_relative_gibbs_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::nondimensional_relative_gibbs_free_energy_per_link(&self.number_of_links, &nondimensional_potential_distance, &nondimensional_potential_stiffness)).into_pyarray(py)
    }
}