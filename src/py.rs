use pyo3::prelude::*;

/// A Python module implemented in Rust.
#[pymodule]
pub fn polymers(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    super::physics::py::register_module(py, m)?;
    Ok(())
}
