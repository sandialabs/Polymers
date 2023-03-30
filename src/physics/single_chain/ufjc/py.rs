use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let ufjc = PyModule::new(py, "ufjc")?;
    super::lennard_jones::py::register_module(py, ufjc)?;
    super::log_squared::py::register_module(py, ufjc)?;
    super::morse::py::register_module(py, ufjc)?;
    parent_module.add_submodule(ufjc)?;
    Ok(())
}
