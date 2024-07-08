use pyo3::prelude::*;
use numpy::
{
    PyArray,
    PyArrayDyn,
    PyReadonlyArrayDyn
};
use super::
{
    NUMGRID,
    init
};

pub fn register_module(_py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    parent_module.add_class::<BucheSilberstein>()?;
    Ok(())
}

/// The Buche-Silberstein hyperelastic damage constitutive model.
///
/// #. Helmholtz method for both the Helmholtz free energy and the equilibrium distribution.
/// #. Gibbs-Legendre method for both the Helmholtz free energy and the equilibrium distribution.
/// #. Gibbs-Legendre for the Helmholtz free energy and a Gaussian equilibrium distribution.
#[pyclass]
pub struct BucheSilberstein
{
    element: [[f64; NUMGRID]; NUMGRID],
    factor: f64,
    grid: [f64; NUMGRID],
    method: u8,
    normalization: f64,

    /// The nondimensional critical extension which irreversibly breaks chains.
    #[pyo3(get)]
    pub nondimensional_critical_extension: f64,

    /// The nondimensional stiffness of each link in a chain.
    #[pyo3(get)]
    pub nondimensional_link_stiffness: f64,

    /// The number of links in a chain.
    #[pyo3(get)]
    pub number_of_links: u8,

    /// The volumetric swelling ratio.
    #[pyo3(get)]
    pub swelling_ratio: f64
}

#[pymethods]
impl BucheSilberstein
{
    #[new]
    pub fn init(method: u8, nondimensional_critical_extension: f64, nondimensional_link_stiffness: f64, number_of_links: u8, swelling_ratio: f64) -> Self
    {
        let (element, factor, grid, normalization) = init(method, nondimensional_critical_extension, nondimensional_link_stiffness, number_of_links, swelling_ratio);
        BucheSilberstein
        {
            element,
            factor,
            grid,
            method,
            normalization,
            nondimensional_critical_extension,
            nondimensional_link_stiffness,
            number_of_links,
            swelling_ratio
        }
    }
    /// The nondimensional Cauchy stress as a function of stretch in uniaxial tension,
    ///
    /// .. math::
    ///     \beta\sigma_{11}(t)/n = 2\pi\int_0^\infty r\,dr\int_0^\infty dz\,P^\mathrm{eq}(\gamma_0)\Theta(\gamma_0; t)\,\frac{\eta(\gamma)}{\gamma}\,\left(2z^2-r^2\right),
    ///
    /// where :math:`\gamma=\sqrt{z^2+r^2}` and :math:`\gamma_0(t)=\sqrt{z^2/F_{11}^2(t)+r^2F_{11}(t)}`. The chain damage function is
    ///
    /// .. math ::
    ///     \Theta(\gamma_0; t) = \begin{cases}
    ///         1, & {}_{(t)}\gamma(s) \leq \gamma_\mathrm{c}~~\forall s\in[0, t], \\
    ///         0, & \mathrm{otherwise},
    ///     \end{cases}
    ///
    /// where :math:`{}_{(t)}\gamma(s)=\sqrt{z^2F_{11}^2(s)/F_{11}^2(t)+r^2F_{11}(t)/F_{11}(s)}`.
    ///
    /// Args:
    ///     stretch (numpy.ndarray): The applied stretch history :math:`F_{11}(t)`.
    /// 
    /// Returns:
    ///     tuple:
    ///       - (*numpy.ndarray*) -
    ///         The nondimensional Cauchy stress :math:`\beta\sigma_{11}(t)/n`.
    ///       - (*numpy.ndarray*) -
    ///         The total probability of intact chains :math:`P^\mathrm{tot}(t)`.
    ///
    pub fn uniaxial_tension<'py>(&self, py: Python<'py>, stretches: PyReadonlyArrayDyn<f64>) -> (&'py PyArrayDyn<f64>, &'py PyArrayDyn<f64>)
    {
        let mut maximum_previous_stretch = 1.0;
        let results = stretches.as_array().iter().map(|stretch|{
            maximum_previous_stretch =
            if maximum_previous_stretch < *stretch {
                *stretch
            } else {
                maximum_previous_stretch
            };
            super::uniaxial_tension(
                &self.element, &self.factor, &self.grid, &self.normalization, &self.method,
                &self.nondimensional_critical_extension, &self.nondimensional_link_stiffness,
                &self.number_of_links, &self.swelling_ratio,
                stretch, &maximum_previous_stretch
            )
        }).collect::<Vec<[f64; 2]>>();
        let results_1 = results.iter().map(|result| result[0]).collect();
        let results_2 = results.iter().map(|result| result[1]).collect();
        (PyArray::from_vec(py, results_1).to_dyn(), PyArray::from_vec(py, results_2).to_dyn())
    }
}
