use pyo3::prelude::*;
use numpy::
{
    IntoPyArray,
    PyArrayDyn,
    PyReadonlyArrayDyn
};

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let legendre = PyModule::new(py, "legendre")?;
    parent_module.add_submodule(legendre)?;
    legendre.add_class::<WLC>()?;
    Ok(())
}

/// The worm-like chain (WLC) model thermodynamics in the isometric ensemble approximated using a Legendre transformation.
#[pyclass]
#[derive(Copy, Clone)]
pub struct WLC
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

    /// The persistance length of the chain in units of nm.
    #[pyo3(get)]
    pub persistance_length: f64,

    nondimensional_persistance_length: f64
}

#[pymethods]
impl WLC
{
    #[new]
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, persistance_length: f64) -> Self
    {
        WLC
        {
            hinge_mass,
            link_length,
            number_of_links,
            persistance_length,
            nondimensional_persistance_length: persistance_length/(number_of_links as f64)/link_length
        }
    }
    /// The Gibbs free energy as a function of the applied end-to-end length and temperature,
    ///
    /// .. math::
    ///     \varphi(\xi, T) \sim \psi(\xi, T) - \xi f(\xi, T) \quad \text{for } N_b\gg 1,
    ///
    /// where :math:`f(\xi, T)` is given by the Legendre transformation approximation above.
    ///
    /// Args:
    ///     end_to_end_length (numpy.ndarray): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The Gibbs free energy :math:`\varphi`.
    ///
    pub fn gibbs_free_energy<'py>(&self, py: Python<'py>, end_to_end_length: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.persistance_length, &end_to_end_length, &temperature)).into_pyarray(py)
    }
    /// The Gibbs free energy per link as a function of the applied end-to-end length and temperature.
    ///
    /// Args:
    ///     end_to_end_length (numpy.ndarray): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The Gibbs free energy per link :math:`\varphi/N_b`.
    ///
    pub fn gibbs_free_energy_per_link<'py>(&self, py: Python<'py>, end_to_end_length: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.persistance_length, &end_to_end_length, &temperature)).into_pyarray(py)
    }
    /// The relative Gibbs free energy as a function of the applied end-to-end length and temperature.
    ///
    /// Args:
    ///     end_to_end_length (numpy.ndarray): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The relative Gibbs free energy :math:`\Delta\varphi\equiv\varphi(\xi,T)-\varphi(0,T)`.
    ///
    pub fn relative_gibbs_free_energy<'py>(&self, py: Python<'py>, end_to_end_length: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::relative_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.persistance_length, &end_to_end_length, &temperature)).into_pyarray(py)
    }
    /// The relative Gibbs free energy per link as a function of the applied end-to-end length and temperature.
    ///
    /// Args:
    ///     end_to_end_length (numpy.ndarray): The end-to-end length :math:`\xi`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The relative Gibbs free energy per link :math:`\Delta\varphi/N_b`.
    ///
    pub fn relative_gibbs_free_energy_per_link<'py>(&self, py: Python<'py>, end_to_end_length: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        end_to_end_length.as_array().mapv(|end_to_end_length: f64| super::relative_gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.persistance_length, &end_to_end_length, &temperature)).into_pyarray(py)
    }
    /// The nondimensional Gibbs free energy as a function of the applied nondimensional end-to-end length per link and temperature.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (numpy.ndarray): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional Gibbs free energy :math:`N_b\varrho=\beta\varphi`.
    ///
    pub fn nondimensional_gibbs_free_energy<'py>(&self, py: Python<'py>, nondimensional_end_to_end_length_per_link: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::nondimensional_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.nondimensional_persistance_length, &nondimensional_end_to_end_length_per_link, &temperature)).into_pyarray(py)
    }
    /// The nondimensional Gibbs free energy per link as a function of the applied nondimensional end-to-end length per link and temperature.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (numpy.ndarray): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional Gibbs free energy per link :math:`\varrho\equiv\beta\varphi/N_b`.
    ///
    pub fn nondimensional_gibbs_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_end_to_end_length_per_link: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::nondimensional_gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.nondimensional_persistance_length, &nondimensional_end_to_end_length_per_link, &temperature)).into_pyarray(py)
    }
    /// The nondimensional relative Gibbs free energy as a function of the applied nondimensional end-to-end length per link.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (numpy.ndarray): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Gibbs free energy :math:`\beta\Delta\varphi=N_b\Delta\varrho`.
    ///
    pub fn nondimensional_relative_gibbs_free_energy<'py>(&self, py: Python<'py>, nondimensional_end_to_end_length_per_link: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::nondimensional_relative_gibbs_free_energy(&self.number_of_links, &self.nondimensional_persistance_length, &nondimensional_end_to_end_length_per_link)).into_pyarray(py)
    }
    /// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional end-to-end length per link.
    ///
    /// Args:
    ///     nondimensional_end_to_end_length_per_link (numpy.ndarray): The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Gibbs free energy per link :math:`\Delta\varrho\equiv\beta\Delta\varphi/N_b`.
    ///
    pub fn nondimensional_relative_gibbs_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_end_to_end_length_per_link: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        nondimensional_end_to_end_length_per_link.as_array().mapv(|nondimensional_end_to_end_length_per_link: f64| super::nondimensional_relative_gibbs_free_energy_per_link(&self.number_of_links, &self.nondimensional_persistance_length, &nondimensional_end_to_end_length_per_link)).into_pyarray(py)
    }
}
