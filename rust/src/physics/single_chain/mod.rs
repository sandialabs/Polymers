use crate::math::Math;

pub mod test;
pub mod test_macros;
pub mod fjc;

pub trait Isometric
{
    fn init(number_of_links: u16, link_length: f64) -> Self;
    fn force<T>(&self, end_to_end_length: &T) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>;
}

pub trait Isotensional
{
    fn init(number_of_links: u16, link_length: f64) -> Self;
    fn end_to_end_length<T>(&self, force: &T, inverse_temperature: f64) -> T;
    fn end_to_end_length_per_link<T>(&self, force: &T, inverse_temperature: f64) -> T;
    fn nondimensional_end_to_end_length<T>(&self, nondimensional_force: &T) -> T;
    fn nondimensional_end_to_end_length_per_link<T>(&self, nondimensional_force: &T) -> T;
    fn relative_gibbs_free_energy<T>(&self, force: &T, inverse_temperature: f64) -> T;
    fn relative_gibbs_free_energy_per_link<T>(&self, force: &T, inverse_temperature: f64) -> T;
    // fn nondimensional_gibbs_free_energy<T>(&self, nondimensional_force: &T) -> T;
    // fn nondimensional_gibbs_free_energy_per_link<T>(&self, nondimensional_force: &T) -> T;
    fn nondimensional_relative_gibbs_free_energy<T>(&self, nondimensional_force: &T) -> T;
    fn nondimensional_relative_gibbs_free_energy_per_link<T>(&self, nondimensional_force: &T) -> T;
}