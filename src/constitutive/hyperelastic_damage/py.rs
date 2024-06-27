use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let hyperelastic_damage = PyModule::new(py, "hyperelastic_damage")?;
    super::buche_silberstein::py::register_module(py, hyperelastic_damage)?;
    parent_module.add_submodule(hyperelastic_damage)?;
    Ok(())
}