use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let isometric = PyModule::new(py, "isometric")?;
    parent_module.add_submodule(isometric)?;
    isometric.add_class::<FJC>()?;
    Ok(())
}

#[pyclass]
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub methods: super::FJC
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
            methods: super::FJC::init(number_of_links, link_length, hinge_mass)
        }
    }
    fn force(&self, end_to_end_length_per_link: f64, temperature: f64) -> PyResult<f64>
    {
        Ok(self.methods.force(&end_to_end_length_per_link, &temperature))
    }
    fn nondimensional_force(&self, nondimensional_end_to_end_length_per_link: f64) -> PyResult<f64>
    {
        Ok(self.methods.nondimensional_force(&nondimensional_end_to_end_length_per_link))
    }
}