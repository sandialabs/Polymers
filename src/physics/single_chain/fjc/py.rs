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
/// The structure of the FJC model.
pub struct FJC
{
    #[pyo3(get)]
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,
    
    #[pyo3(get)]
    /// The length of each link in the chain in units of nm.
    pub link_length: f64,
    
    #[pyo3(get)]
    /// The number of links in the chain.
    pub number_of_links: u8,
    
    #[pyo3(get)]
    /// The thermodynamic functions of the model.
    pub thermodynamics: super::thermodynamics::py::FJC
}

#[pymethods]
/// The implemented functionality of the FJC model.
impl FJC
{
    #[new]
    /// Initializes and returns an instance of the FJC model.
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
    fn __call__(&self) -> PyResult<(u8, f64, f64)>
    {
        Ok((self.number_of_links, self.link_length, self.hinge_mass))
    }
    fn __repr__(&self) -> PyResult<String>
    {
        Ok(format!("FJC(number_of_links={N}, link_length={l}, hinge_mass={m})", N=self.number_of_links, l=self.link_length, m=self.hinge_mass).into())
    }
}