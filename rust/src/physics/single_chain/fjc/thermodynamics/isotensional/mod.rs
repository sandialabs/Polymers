use crate::math::
{
    Math,
    langevin,
    ln_sinhc
};

pub mod test;

pub struct Isotensional
{
    pub number_of_links: u16,
    pub number_of_links_f64: f64,
    pub link_length: f64,
    pub legendre: Legendre
}

impl Isotensional
{
    pub fn init(number_of_links: u16, link_length: f64) -> Isotensional
    {
        Isotensional
        {
            number_of_links: number_of_links,
            number_of_links_f64: number_of_links as f64,
            link_length: link_length,
            legendre: Legendre::init(number_of_links, link_length)
        }
    }
    pub fn end_to_end_length<T>(&self, force: &T, inverse_temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Mul<f64, Output = T>
    {
        self.end_to_end_length_per_link(force, inverse_temperature)*self.number_of_links_f64
    }
    pub fn end_to_end_length_per_link<T>(&self, force: &T, inverse_temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Mul<f64, Output = T>
    {
        let nondimensional_force = &(*force*inverse_temperature*self.link_length);
        self.nondimensional_end_to_end_length_per_link(nondimensional_force)*self.link_length
    }
    pub fn nondimensional_end_to_end_length<T>(&self, nondimensional_force: &T) -> T
    where T:
        Math<T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_end_to_end_length_per_link(nondimensional_force)*self.number_of_links_f64
    }
    pub fn nondimensional_end_to_end_length_per_link<T>(&self, nondimensional_force: &T) -> T
    where T:
        Math<T>
    {
        langevin(nondimensional_force)
    }
    // pub fn gibbs_free_energy<T>(&self, force: &T, inverse_temperature: f64)
    // {
        
    // }
    // pub fn gibbs_free_energy_per_link<T>(&self, force: &T, inverse_temperature: f64)
    // {

    // }
    pub fn relative_gibbs_free_energy<T>(&self, force: &T, inverse_temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.relative_gibbs_free_energy_per_link(force, inverse_temperature)*self.number_of_links_f64
    }
    pub fn relative_gibbs_free_energy_per_link<T>(&self, force: &T, inverse_temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        let nondimensional_force = &(*force*inverse_temperature*self.link_length);
        self.nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force)/inverse_temperature
    }
    // pub fn nondimensional_gibbs_free_energy<T>(&self, nondimensional_force: &T)
    // {

    // }
    // pub fn nondimensional_gibbs_free_energy_per_link<T>(&self, nondimensional_force: &T)
    // {

    // }
    pub fn nondimensional_relative_gibbs_free_energy<T>(&self, nondimensional_force: &T) -> T
    where T:
        Math<T> +
        std::ops::Neg<Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force)*self.number_of_links_f64
    }
    pub fn nondimensional_relative_gibbs_free_energy_per_link<T>(&self, nondimensional_force: &T) -> T
    where T:
        Math<T> +
        std::ops::Neg<Output = T>
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
            number_of_links: number_of_links,
            number_of_links_f64: number_of_links as f64,
            link_length: link_length
        }
    }
    // pub fn helmholtz_free_energy<T>(&self, nondimensional_force: &T)
    // {

    // }
    // pub fn helmholtz_free_energy_per_link<T>(&self, nondimensional_force: &T)
    // {

    // }
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
        self.nondimensional_relative_helmholtz_free_energy(nondimensional_force)/inverse_temperature
    }
    // pub fn nondimensional_helmholtz_free_energy<T>(&self, nondimensional_force: &T)
    // {

    // }
    // pub fn nondimensional_helmholtz_free_energy_per_link<T>(&self, nondimensional_force: &T)
    // {

    // }
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