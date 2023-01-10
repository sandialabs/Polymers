use pyo3::prelude::*;

/// Python API for the library.

#[pymodule]
pub fn polymers(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    super::physics::py::register_module(py, m)?;
    Ok(())
}
