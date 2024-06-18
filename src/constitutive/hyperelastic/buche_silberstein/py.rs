use pyo3::prelude::*;
use numpy::
{
    IntoPyArray,
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

/// The Buche-Silberstein hyperelastic constitutive model.
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

    /// The nondimensional stiffness of each link in a chain.
    #[pyo3(get)]
    pub nondimensional_link_stiffness: f64,

    /// The number of links in a chain.
    #[pyo3(get)]
    pub number_of_links: u8
}

#[pymethods]
impl BucheSilberstein
{
    #[new]
    pub fn init(method: u8, nondimensional_link_stiffness: f64, number_of_links: u8) -> Self
    {
        init(method, nondimensional_link_stiffness, number_of_links);
        let (element, factor, grid, normalization) = init(method, nondimensional_link_stiffness, number_of_links);
        BucheSilberstein
        {
            element,
            factor,
            grid,
            method,
            normalization,
            nondimensional_link_stiffness,
            number_of_links
        }
    }
    /// The nondimensional Cauchy stress as a function of stretch in uniaxial tension,
    ///
    /// .. math::
    ///     \beta\sigma_{11}/n = 2\pi\int_0^\infty r\,dr\int_0^\infty dz\,P^\mathrm{eq}(\gamma_0)\,\frac{\eta(\gamma)}{\gamma}\,\left(2z^2-r^2\right),
    ///
    /// where :math:`\gamma=\sqrt{z^2+r^2}` and :math:`\gamma_0=\sqrt{\left(z/F_{11}\right)^2+F_{11}r^2}`.
    ///
    /// Args:
    ///     stretch (numpy.ndarray): The applied stretch :math:`F_{11}`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional Cauchy stress :math:`\beta\sigma_{11}/n`.
    ///
    pub fn uniaxial_tension<'py>(&self, py: Python<'py>, stretch: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        stretch.as_array().mapv(|stretch: f64|
            super::uniaxial_tension(
                &self.element, &self.factor, &self.grid, &self.normalization, &self.method,
                &self.nondimensional_link_stiffness, &self.number_of_links, &stretch
            )
        ).into_pyarray(py)
    }
    /// The nondimensional Cauchy stress as a function of stretch in equibiaxial tension,
    ///
    /// .. math::
    ///     \beta\sigma_{11}/n = 2\pi\int_0^\infty r\,dr\int_0^\infty dz\,P^\mathrm{eq}(\gamma_0)\,\frac{\eta(\gamma)}{\gamma}\,\left(r^2-2z^2\right),
    ///
    /// where :math:`\gamma=\sqrt{z^2+r^2}` and :math:`\gamma_0=\sqrt{\left(F_{11}^2z\right)^2+\left(r/F_{11}\right)^2}`.
    ///
    /// Args:
    ///     stretch (numpy.ndarray): The applied stretch :math:`F_{11}`.
    /// 
    /// Returns:
    ///     numpy.ndarray: The nondimensional Cauchy stress :math:`\beta\sigma_{11}/n`.
    ///
    pub fn equibiaxial_tension<'py>(&self, py: Python<'py>, stretch: PyReadonlyArrayDyn<f64>) -> &'py PyArrayDyn<f64>
    {
        stretch.as_array().mapv(|stretch: f64|
            super::equibiaxial_tension(
                &self.element, &self.factor, &self.grid, &self.normalization, &self.method,
                &self.nondimensional_link_stiffness, &self.number_of_links, &stretch
            )
        ).into_pyarray(py)
    }
}