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
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub thermodynamics: super::thermodynamics::py::FJC
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
            thermodynamics: super::thermodynamics::py::FJC::init(number_of_links, link_length, hinge_mass)
        }
    }
    fn __call__(&self) -> PyResult<(u8, f64, f64)>
    {
        Ok((self.number_of_links, self.link_length, self.hinge_mass))
    }
    fn __repr__(&self) -> PyResult<String>
    {
        Ok(format!("FJC(number_of_links={N}, link_length={l}, hinge_mass={m})", N=self.number_of_links, l=self.link_length, m=self.hinge_mass).into())
    }
}