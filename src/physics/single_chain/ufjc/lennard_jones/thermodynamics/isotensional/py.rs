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
    let isotensional = PyModule::new(py, "isotensional")?;
    super::asymptotic::py::register_module(py, isotensional)?;
    super::legendre::py::register_module(py, isotensional)?;
    parent_module.add_submodule(isotensional)?;
    isotensional.add_class::<LENNARDJONESFJC>()?;
    Ok(())
}

/// The Lennard-Jones link potential freely-jointed chain (log-squared-FJC) model thermodynamics in the isotensional ensemble.
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
    pub link_stiffness: f64,

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using an asymptotic approach.
    #[pyo3(get)]
    pub asymptotic: super::asymptotic::py::LENNARDJONESFJC,

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using a Legendre transformation.
    #[pyo3(get)]
    pub legendre: super::legendre::py::LENNARDJONESFJC
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
            link_stiffness,
            asymptotic: super::asymptotic::py::LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness),
            legendre: super::legendre::py::LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness)
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
        force.as_array().mapv(|force: f64| super::end_to_end_length(&self.number_of_links, &self.link_length, &self.link_stiffness, &force, &temperature)).into_pyarray(py)
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
        force.as_array().mapv(|force: f64| super::end_to_end_length_per_link(&self.link_length, &self.link_stiffness, &force, &temperature)).into_pyarray(py)
    }
    /// The expected nondimensional end-to-end length as a function of the applied nondimensional force.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional end-to-end length :math:`N_b\gamma=\xi/\ell_b`.
    ///
    pub fn nondimensional_end_to_end_length<'py>(&self, py: Python<'py>, nondimensional_force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_end_to_end_length(&self.number_of_links, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_force)).into_pyarray(py)
    }
    /// The expected nondimensional end-to-end length per link as a function of the applied nondimensional force, given by :footcite:t:`buche2022freely` as
    ///
    /// .. math::
    ///     \gamma(\eta) = -\frac{\partial}{\partial\eta}\,\ln\left[\int \frac{\sinh(s\eta)}{s\eta}\,e^{-\beta u(s)}s^2\,ds\right].
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional end-to-end length per link :math:`\gamma\equiv \xi/N_b\ell_b`.
    ///
    pub fn nondimensional_end_to_end_length_per_link<'py>(&self, py: Python<'py>, nondimensional_force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_end_to_end_length_per_link(&(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_force)).into_pyarray(py)
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
        force.as_array().mapv(|force: f64| super::gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.link_stiffness, &force, &temperature)).into_pyarray(py)
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
        force.as_array().mapv(|force: f64| super::gibbs_free_energy_per_link(&self.link_length, &self.hinge_mass, &self.link_stiffness, &force, &temperature)).into_pyarray(py)
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
        force.as_array().mapv(|force: f64| super::relative_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.link_stiffness, &force, &temperature)).into_pyarray(py)
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
        force.as_array().mapv(|force: f64| super::relative_gibbs_free_energy_per_link(&self.link_length, &self.link_stiffness, &force, &temperature)).into_pyarray(py)
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
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_force, &temperature)).into_pyarray(py)
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
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_gibbs_free_energy_per_link(&self.link_length, &self.hinge_mass, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_force, &temperature)).into_pyarray(py)
    }
    /// The nondimensional relative Gibbs free energy as a function of the applied nondimensional force.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Gibbs free energy :math:`\beta\Delta\varphi=N_b\Delta\varrho`.
    ///
    pub fn nondimensional_relative_gibbs_free_energy<'py>(&self, py: Python<'py>, nondimensional_force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_relative_gibbs_free_energy(&self.number_of_links, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_force)).into_pyarray(py)
    }
    /// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional force, given by :footcite:t:`buche2022freely` as
    ///
    /// .. math::
    ///     \Delta\varrho(\eta) = \ln\left[\int \frac{\sinh(s\eta)}{s\eta}\,e^{-\beta u(s)}s^2\,ds\right]
    ///                         - \ln\left[\int e^{-\beta u(s)}s^2\,ds\right],
    ///
    /// where the nondimensional link potential :math:`\beta u` is given by :footcite:t:`jones1924determinationii` as
    ///
    /// .. math::
    ///     \beta u(\lambda) = \varepsilon\left[\frac{1}{\lambda^{12}} - \frac{2}{\lambda^6}\right],
    ///
    /// where :math:`\varepsilon\equiv\beta u_b=\kappa/72` is the nondimensional potential energy scale,
    /// :math:`\kappa\equiv\beta k_b\ell_b^2` is the nondimensional link stiffness,
    /// and :math:`\lambda\equiv\ell/\ell_b` is the nondimensional link stretch.
    ///
    /// Args:
    ///     nondimensional_force (numpy.ndarray): The nondimensional force :math:`\eta\equiv\beta f\ell_b`.
    ///     temperature (float): The temperature :math:`T`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional relative Gibbs free energy per link :math:`\Delta\varrho\equiv\beta\Delta\varphi/N_b`.
    ///
    pub fn nondimensional_relative_gibbs_free_energy_per_link<'py>(&self, py: Python<'py>, nondimensional_force: PyReadonlyArrayDyn<f64>, temperature: f64) -> &'py PyArrayDyn<f64>
    {
        nondimensional_force.as_array().mapv(|nondimensional_force: f64| super::nondimensional_relative_gibbs_free_energy_per_link(&(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &nondimensional_force)).into_pyarray(py)
    }
}