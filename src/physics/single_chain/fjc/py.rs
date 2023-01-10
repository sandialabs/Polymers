use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let fjc = PyModule::new(py, "fjc")?;
    super::thermodynamics::py::register_module(py, fjc)?;
    parent_module.add_submodule(fjc)?;
    fjc.add_class::<FJC>()?;
    Ok(())
}

#[pyclass]
pub struct FJC
{
    #[pyo3(get)]
    pub hinge_mass: f64,
    
    #[pyo3(get)]
    pub link_length: f64,
    
    #[pyo3(get)]
    pub number_of_links: u8,
    
    #[pyo3(get)]
    pub thermodynamics: super::thermodynamics::py::FJC
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
            thermodynamics: super::thermodynamics::py::FJC::init(number_of_links, link_length, hinge_mass)
        }
    }
}
