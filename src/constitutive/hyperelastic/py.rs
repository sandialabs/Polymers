use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let hyperelastic = PyModule::new(py, "hyperelastic")?;
    super::buche_silberstein::py::register_module(py, hyperelastic)?;
    parent_module.add_submodule(hyperelastic)?;
    Ok(())
}