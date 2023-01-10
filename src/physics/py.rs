//!
//! asdfasdfasdfs
use pyo3::prelude::*;

/// A Python module implemented in Rust (1).
#[pymodule]
/// A Python module implemented in Rust (2).
pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let physics = PyModule::new(py, "physics")?;
    super::single_chain::py::register_module(py, physics)?;
    parent_module.add_submodule(physics)?;
    Ok(())
}
