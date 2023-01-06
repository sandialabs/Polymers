use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let single_chain = PyModule::new(py, "single_chain")?;
    super::fjc::py::register_module(py, single_chain)?;
    parent_module.add_submodule(single_chain)?;
    Ok(())
}
