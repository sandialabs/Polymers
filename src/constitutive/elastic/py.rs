use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let elastic = PyModule::new(py, "elastic")?;
    super::neo_hookean::py::register_module(py, elastic)?;
    parent_module.add_submodule(elastic)?;
    Ok(())
}
