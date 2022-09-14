use crate::math::
{
    Math,
    langevin,
    ln_sinhc
};
use crate::physics::single_chain::Isotensional;

pub mod test;

pub struct FJC
{
    pub link_length: f64,
    pub number_of_links: u16,
    pub number_of_links_f64: f64,
    pub legendre: Legendre
}

impl Isotensional for FJC
{
    fn init(number_of_links: u16, link_length: f64) -> FJC
    {
        FJC
        {
            link_length: link_length,
            number_of_links: number_of_links,
            number_of_links_f64: number_of_links as f64,
            legendre: Legendre::init(number_of_links, link_length)
        }
    }
    fn end_to_end_length<T>(&self, force: &T, inverse_temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.end_to_end_length_per_link(force, inverse_temperature)*self.number_of_links_f64
    }
    fn end_to_end_length_per_link<T>(&self, force: &T, inverse_temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        let nondimensional_force = &(*force*inverse_temperature*self.link_length);
        self.nondimensional_end_to_end_length_per_link(nondimensional_force)*self.link_length
    }
    fn nondimensional_end_to_end_length<T>(&self, nondimensional_force: &T) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_end_to_end_length_per_link(nondimensional_force)*self.number_of_links_f64
    }
    fn nondimensional_end_to_end_length_per_link<T>(&self, nondimensional_force: &T) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        langevin(nondimensional_force)
    }
    // pub fn gibbs_free_energy<T>(&self, force: &T, inverse_temperature: f64)
        // where T:
        // Math<T> +
        // std::marker::Copy +
        // std::ops::Neg<Output = T> +
        // std::ops::Mul<T, Output = T> +
        // std::ops::Sub<T, Output = T> +
        // std::ops::Div<f64, Output = T> +
        // std::ops::Mul<f64, Output = T>
    // {}
    // pub fn gibbs_free_energy_per_link<T>(&self, force: &T, inverse_temperature: f64)
        // where T:
        // Math<T> +
        // std::marker::Copy +
        // std::ops::Neg<Output = T> +
        // std::ops::Mul<T, Output = T> +
        // std::ops::Sub<T, Output = T> +
        // std::ops::Div<f64, Output = T> +
        // std::ops::Mul<f64, Output = T>
    // {}
    fn relative_gibbs_free_energy<T>(&self, force: &T, inverse_temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.relative_gibbs_free_energy_per_link(force, inverse_temperature)*self.number_of_links_f64
    }
    fn relative_gibbs_free_energy_per_link<T>(&self, force: &T, inverse_temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        let nondimensional_force = &(*force*inverse_temperature*self.link_length);
        self.nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force)/inverse_temperature
    }
    // fn nondimensional_gibbs_free_energy<T>(&self, nondimensional_force: &T) -> T
        // where T:
        // Math<T> +
        // std::marker::Copy +
        // std::ops::Neg<Output = T> +
        // std::ops::Mul<T, Output = T> +
        // std::ops::Sub<T, Output = T> +
        // std::ops::Div<f64, Output = T> +
        // std::ops::Mul<f64, Output = T>
    // {}
    // fn nondimensional_gibbs_free_energy_per_link<T>(&self, nondimensional_force: &T) -> T
        // where T:
        // Math<T> +
        // std::marker::Copy +
        // std::ops::Neg<Output = T> +
        // std::ops::Mul<T, Output = T> +
        // std::ops::Sub<T, Output = T> +
        // std::ops::Div<f64, Output = T> +
        // std::ops::Mul<f64, Output = T>
    // {}
    fn nondimensional_relative_gibbs_free_energy<T>(&self, nondimensional_force: &T) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force)*self.number_of_links_f64
    }
    fn nondimensional_relative_gibbs_free_energy_per_link<T>(&self, nondimensional_force: &T) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        -ln_sinhc(nondimensional_force)
    }
}

pub struct Legendre
{
    pub number_of_links: u16,
    pub number_of_links_f64: f64,
    pub link_length: f64
}

impl Legendre
{
    pub fn init(number_of_links: u16, link_length: f64) -> Legendre
    {
        Legendre
        {
            link_length: link_length,
            number_of_links: number_of_links,
            number_of_links_f64: number_of_links as f64
        }
    }
    // pub fn helmholtz_free_energy<T>(&self, nondimensional_force: &T)
    // {}
    // pub fn helmholtz_free_energy_per_link<T>(&self, nondimensional_force: &T)
    // {}
    pub fn relative_helmholtz_free_energy<T>(&self, force: &T, inverse_temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.relative_helmholtz_free_energy_per_link(force, inverse_temperature)*self.number_of_links_f64
    }
    pub fn relative_helmholtz_free_energy_per_link<T>(&self, force: &T, inverse_temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        let nondimensional_force = &(*force*inverse_temperature*self.link_length);
        self.nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_force)/inverse_temperature
    }
    // pub fn nondimensional_helmholtz_free_energy<T>(&self, nondimensional_force: &T)
    // {}
    // pub fn nondimensional_helmholtz_free_energy_per_link<T>(&self, nondimensional_force: &T)
    // {}
    pub fn nondimensional_relative_helmholtz_free_energy<T>(&self, nondimensional_force: &T) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_force)*self.number_of_links_f64
    }
    pub fn nondimensional_relative_helmholtz_free_energy_per_link<T>(&self, nondimensional_force: &T) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T>
    {
        *nondimensional_force*langevin(nondimensional_force) - ln_sinhc(nondimensional_force)
    }
}