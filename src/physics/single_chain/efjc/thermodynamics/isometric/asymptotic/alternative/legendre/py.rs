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
    legendre.add_class::<EFJC>()?;
    Ok(())
}

/// The extensible freely-jointed chain (EFJC) model thermodynamics in the isometric ensemble approximated using an alternative asymptotic approach and a Legendre transformation.
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
    pub link_stiffness: f64
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
            link_stiffness
        }
    }
    /// The expected force as a function of the applied end-to-end length and temperature.
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
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::force(&self.number_of_links, &self.link_length, &self.link_stiffness, &end_to_end_length, &temperature)).into_pyarray(py)
    }
    /// The expected nondimensional force as a function of the applied nondimensional end-to-end length per link.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (numpy.ndarray): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///
    pub fn nondimensional_force<'py>(&self, py: Python<'py>, nondimensional_end_to_end_length_per_link: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::nondimensional_force(&(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_end_to_end_length_per_link)).into_pyarray(py)
    }
    /// The Helmholtz free energy as a function of the applied end-to-end length and temperature,
    ///
    /// .. math::
    ///     \psi(\xi, T) \sim \varphi\left[f(\xi, T)\right] + \xi f(\xi, T) \quad \text{for } N_b\gg 1,
    ///
    /// where :math:`f(\xi, T)` is given by the Legendre transformation approximation above.
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
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.link_stiffness, &end_to_end_length, &temperature)).into_pyarray(py)
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
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.link_stiffness, &end_to_end_length, &temperature)).into_pyarray(py)
    }
    /// The relative Helmholtz free energy as a function of the applied end-to-end length and temperature.
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
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::relative_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.link_stiffness, &end_to_end_length, &temperature)).into_pyarray(py)
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
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::relative_helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.link_stiffness, &end_to_end_length, &temperature)).into_pyarray(py)
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
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::nondimensional_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_end_to_end_length_per_link, &temperature)).into_pyarray(py)
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
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::nondimensional_helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_end_to_end_length_per_link, &temperature)).into_pyarray(py)
    }
    /// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional end-to-end length per link.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (numpy.ndarray): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Helmholtz free energy :math:`\beta\Delta\psi=N_b\Delta\vartheta`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy<'py>(&self, py: Python<'py>, nondimensional_end_to_end_length_per_link: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::nondimensional_relative_helmholtz_free_energy(&self.number_of_links, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_end_to_end_length_per_link)).into_pyarray(py)
    }
    /// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional end-to-end length per link.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (numpy.ndarray): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Helmholtz free energy per link :math:`\Delta\vartheta\equiv\beta\Delta\psi/N_b`.
    ///
    pub fn nondimensional_relative_helmholtz_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_end_to_end_length_per_link: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::nondimensional_relative_helmholtz_free_energy_per_link(&(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_end_to_end_length_per_link)).into_pyarray(py)
    }
}