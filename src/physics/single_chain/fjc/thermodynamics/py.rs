use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let thermodynamics = PyModule::new(py, "thermodynamics")?;
    super::isometric::py::register_module(py, thermodynamics)?;
    parent_module.add_submodule(thermodynamics)?;
    thermodynamics.add_class::<FJC>()?;
    Ok(())
}

#[pyclass]
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub isometric: super::isometric::py::FJC
}

#[pymethods]
impl FJC
{
    #[new]
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> FJC
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            isometric: super::isometric::py::FJC::init(number_of_links, link_length, hinge_mass)
        }
    }
}