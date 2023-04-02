use pyo3::prelude::*;

pub fn register_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()>
{
    let math = PyModule::new(py, "math")?;
    parent_module.add_submodule(math)?;
    math.add_function(wrap_pyfunction!(integrate_1d, math)?)?;
    Ok(())
}

#[pyfunction]
pub fn integrate_1d(f: &dyn Fn(&f64) -> f64, x_min: f64, x_max: f64, num_points: u128) -> f64
{
    super::integrate_1d(f, &x_min, &x_max, &num_points)
}