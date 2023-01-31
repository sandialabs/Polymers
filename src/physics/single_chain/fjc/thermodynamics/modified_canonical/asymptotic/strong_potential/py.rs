use pyo3::prelude::*;
use numpy::
{
    IntoPyArray,
    PyArrayDyn,
    PyReadonlyArrayDyn
};

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
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The force :math:`f`.
    ///
    pub fn force<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).force(&potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
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
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness)).into_pyarray(py)
    }
    /// The helmholtz free energy as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The helmholtz free energy :math:`\psi`.
    ///
    pub fn helmholtz_free_energy<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The helmholtz free energy per link as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The helmholtz free energy per link :math:`\psi/N_b`.
    ///
    pub fn helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The relative helmholtz free energy as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The relative helmholtz free energy :math:`\Delta\psi`.
    ///
    pub fn relative_helmholtz_free_energy<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).relative_helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The relative helmholtz free energy per link as a function of the applied potential distance, potential stiffness, and temperature.
    ///
    /// Args:
    ///     potential_distance (numpy.ndarray): The potential distance.
    ///     potential_stiffness (float): The potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The relative helmholtz free energy per link :math:`\Delta\psi/N_b`.
    ///
    pub fn relative_helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, potential_distance: PyReadonlyArrayDyn<f64>, potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        potential_distance.as_array().mapv(|potential_distance: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).relative_helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The nondimensional helmholtz free energy as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The nondimensional helmholtz free energy :math:`\beta\psi=N_b\vartheta`.
    ///
    pub fn nondimensional_helmholtz_free_energy<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The nondimensional helmholtz free energy per link as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///     temperature (float): The temperature :math:`T`.
    ///
    /// Returns:
    ///     numpy.ndarray: The nondimensional helmholtz free energy per link :math:`\vartheta\equiv\beta\psi/N_b`.
    ///
    pub fn nondimensional_helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature)).into_pyarray(py)
    }
    /// The nondimensional relative helmholtz free energy as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative helmholtz free energy :math:`\beta\Delta\psi=N_b\Delta\vartheta`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_relative_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness)).into_pyarray(py)
    }
    /// The nondimensional relative helmholtz free energy per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    ///
    /// Args:
    ///     nondimensional_potential_distance (numpy.ndarray): The nondimensional potential distance.
    ///     nondimensional_potential_stiffness (float): The nondimensional potential stiffness.
    ///
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative helmholtz free energy per link :math:`\Delta\vartheta\equiv\beta\Delta\psi/N_b`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_potential_distance: PyReadonlyArrayDyn<f64>, nondimensional_potential_stiffness: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_potential_distance.as_array().mapv(|nondimensional_potential_distance: f64| super::FJC::init(self.number_of_links, self.link_length, self.hinge_mass).nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness)).into_pyarray(py)
    }
}