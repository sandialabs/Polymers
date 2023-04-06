use pyo3::prelude::*;
use numpy::
{
    IntoPyArray,
    PyArrayDyn,
    PyReadonlyArrayDyn
};
use crate::physics::BOLTZMANN_CONSTANT;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let legendre = PyModule::new(py, "legendre")?;
    parent_module.add_submodule(legendre)?;
    legendre.add_class::<LENNARDJONESFJC>()?;
    Ok(())
}

/// The Lennard-Jones link potential freely-jointed chain (Lennard-Jones-FJC) model thermodynamics in the isotensional ensemble approximated using a Legendre transformation.
#[pyclass]
#[derive(Copy, Clone)]
pub struct LENNARDJONESFJC
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
    pub link_stiffness: f64
}

#[pymethods]
impl LENNARDJONESFJC
{
    #[new]
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        LENNARDJONESFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness
        }
    }
    /// The Helmholtz free energy as a function of the applied force and temperature,
    ///
    /// .. math::
    ///     \psi(f, T) \sim \varphi(f, T) + f \xi(f, T) \quad \text{for } N_b\gg 1.
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The Helmholtz free energy :math:`\psi`.
    ///
    pub fn helmholtz_free_energy<'py>(&self, py: Python<'py>, force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        force.as_array().mapv(|force: f64| super::helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.link_stiffness, &force, &temperature)).into_pyarray(py)
    }
    /// The Helmholtz free energy per link as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The Helmholtz free energy per link :math:`\psi/N_b`.
    ///
    pub fn helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        force.as_array().mapv(|force: f64| super::helmholtz_free_energy_per_link(&self.link_length, &self.hinge_mass, &self.link_stiffness, &force, &temperature)).into_pyarray(py)
    }
    /// The relative Helmholtz free energy as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The relative Helmholtz free energy :math:`\Delta\psi\equiv\psi(f,T)-\psi(0,T)`.
    ///
    pub fn relative_helmholtz_free_energy<'py>(&self, py: Python<'py>, force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        force.as_array().mapv(|force: f64| super::relative_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.link_stiffness, &force, &temperature)).into_pyarray(py)
    }
    /// The relative Helmholtz free energy per link as a function of the applied force and temperature.
    ///
    /// Args:
    ///     force (numpy.ndarray): The force :math:`f`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The relative Helmholtz free energy per link :math:`\Delta\psi/N_b`.
    ///
    pub fn relative_helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        force.as_array().mapv(|force: f64| super::relative_helmholtz_free_energy_per_link(&self.link_length, &self.link_stiffness, &force, &temperature)).into_pyarray(py)
    }
    /// The nondimensional Helmholtz free energy as a function of the applied nondimensional force and temperature.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional Helmholtz free energy :math:`\beta\psi=N_b\vartheta`.
    ///
    pub fn nondimensional_helmholtz_free_energy<'py>(&self, py: Python<'py>, nondimensional_force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_force, &temperature)).into_pyarray(py)
    }
    /// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional force and temperature.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional Helmholtz free energy per link :math:`\vartheta\equiv\beta\psi/N_b`.
    ///
    pub fn nondimensional_helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_helmholtz_free_energy_per_link(&self.link_length, &self.hinge_mass, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_force, &temperature)).into_pyarray(py)
    }
    /// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional force.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Helmholtz free energy :math:`\beta\Delta\psi=N_b\Delta\vartheta`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy<'py>(&self, py: Python<'py>, nondimensional_force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_relative_helmholtz_free_energy(&self.number_of_links, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_force)).into_pyarray(py)
    }
    /// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional force.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Helmholtz free energy per link :math:`\Delta\vartheta\equiv\beta\Delta\psi/N_b`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_relative_helmholtz_free_energy_per_link(&(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_force)).into_pyarray(py)
    }
}
