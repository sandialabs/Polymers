use pyo3::prelude::*;

#[pymodule]
pub fn polymers(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    super::physics::py::register_module(py, m)?;
    Ok(())
}
