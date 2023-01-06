use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let physics = PyModule::new(py, "physics")?;
    super::single_chain::py::register_module(py, physics)?;
    parent_module.add_submodule(physics)?;
    Ok(())
}
