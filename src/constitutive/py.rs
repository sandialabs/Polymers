use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let constitutive = PyModule::new(py, "constitutive")?;
    super::hyperelastic::py::register_module(py, constitutive)?;
    parent_module.add_submodule(constitutive)?;
    Ok(())
}
