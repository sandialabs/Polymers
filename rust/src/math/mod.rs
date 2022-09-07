use ndarray::{
    Array,
    Array1,
    Array2,
    Axis,
    stack
};
use ndarray_rand::RandomExt;
use ndarray_rand::rand_distr::Uniform;

pub mod test;

// test more moments and stuff (also finish tests for below)
// implement for a since vector too
// similarly for random_config
// possibly find good optimization between array size and overall samples for Monte Carlo using these
pub fn random_unit_vector(number_of_unit_vectors: u32) -> Array2<f64>
{
    let pi = std::f64::consts::PI;
    let azimuthal_angle: Array1<f64> = 2.0*pi*Array::random(number_of_unit_vectors as usize, Uniform::new(0.0, 1.0));
    let x: Array1<f64> = Array::random(number_of_unit_vectors as usize, Uniform::new(-1.0, 1.0));
    let elevation_angle: Array1<f64> = x.mapv_into(|v| v.acos());
    let azimuthal_sine = azimuthal_angle.to_owned().mapv_into(|v| v.sin());
    let elevation_sine = elevation_angle.to_owned().mapv_into(|v| v.sin());
    let azimuthal_cosine = azimuthal_angle.to_owned().mapv_into(|v| v.cos());
    let elevation_cosine = elevation_angle.to_owned().mapv_into(|v| v.cos());
    stack![Axis(0),azimuthal_cosine*elevation_sine.to_owned(), azimuthal_sine*elevation_sine, elevation_cosine]
}

pub fn sinhc<T>(nondimensional_energy: &T) -> T
where
    T: Math<T>
{
    Math::<T>::sinhc(nondimensional_energy)
}

pub fn ln_sinhc<T>(nondimensional_energy: &T) -> T
where
    T: Math<T>
{
    Math::<T>::ln_sinhc(nondimensional_energy)
}

pub fn langevin<T>(nondimensional_energy: &T) -> T
where
    T: Math<T>
{
    Math::<T>::langevin(nondimensional_energy)
}

pub trait Math<T>
{
    fn sinhc(&self) -> T;
    fn ln_sinhc(&self) -> T;
    fn langevin(&self) -> T;
}

impl Math<f64> for f64
{
    fn sinhc(&self) -> Self
    {
        self.sinh()/self
    }
    fn ln_sinhc(&self) -> Self
    {
        (self.sinh()/self).ln()
    }
    fn langevin(&self) -> Self
    {
        1.0/self.tanh() - 1.0/self
    }
}

impl Math<Array1<f64>> for Array1<f64>
{
    fn sinhc(&self) -> Self
    {
        self.to_owned().mapv_into(|v| v.sinh())/self
    }
    fn ln_sinhc(&self) -> Self
    {
        (self.to_owned().mapv_into(|v| v.sinh())/self).mapv_into(|v| v.ln())
    }
    fn langevin(&self) -> Self
    {
        1.0/self.to_owned().mapv_into(|v| v.tanh()) - 1.0/self
    }
}

impl Math<Array2<f64>> for Array2<f64>
{
    fn sinhc(&self) -> Self
    {
        self.to_owned().mapv_into(|v| v.sinh())/self
    }
    fn ln_sinhc(&self) -> Self
    {
        (self.to_owned().mapv_into(|v| v.sinh())/self).mapv_into(|v| v.ln())
    }
    fn langevin(&self) -> Self
    {
        1.0/self.to_owned().mapv_into(|v| v.tanh()) - 1.0/self
    }
}