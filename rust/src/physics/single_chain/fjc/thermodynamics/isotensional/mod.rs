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
    pub link_length_in_meters: f64,
    pub legendre: Legendre
}

impl Isotensional
{
    pub fn init(number_of_links: u16) -> Isotensional
    {
        Isotensional
        {
            number_of_links: number_of_links,
            number_of_links_f64: number_of_links as f64,
            link_length_in_meters: 1.0,
            legendre: Legendre::init(number_of_links)
        }
    }
    pub fn end_to_end_length_in_meters<T>(&self, force_in_newtons: &T, inverse_temperature_in_inverse_joules: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Mul<f64, Output = T>
    {
        self.end_to_end_length_per_link_in_meters(force_in_newtons, inverse_temperature_in_inverse_joules)*self.number_of_links_f64
    }
    pub fn end_to_end_length_per_link_in_meters<T>(&self, force_in_newtons: &T, inverse_temperature_in_inverse_joules: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Mul<f64, Output = T>
    {
        let nondimensional_force = &(*force_in_newtons*inverse_temperature_in_inverse_joules*self.number_of_links_f64);
        self.nondimensional_end_to_end_length_per_link(nondimensional_force)*self.link_length_in_meters
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
    // pub fn gibbs_free_energy_in_joules<T>(&self, force_in_newtons: &T, inverse_temperature_in_inverse_joules: f64)
    // {
        
    // }
    // pub fn gibbs_free_energy_per_link_in_joules<T>(&self, force_in_newtons: &T, inverse_temperature_in_inverse_joules: f64)
    // {

    // }
    pub fn relative_gibbs_free_energy_in_joules<T>(&self, force_in_newtons: &T, inverse_temperature_in_inverse_joules: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.relative_gibbs_free_energy_per_link_in_joules(force_in_newtons, inverse_temperature_in_inverse_joules)*self.number_of_links_f64
    }
    pub fn relative_gibbs_free_energy_per_link_in_joules<T>(&self, force_in_newtons: &T, inverse_temperature_in_inverse_joules: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        let nondimensional_force = &(*force_in_newtons*inverse_temperature_in_inverse_joules*self.number_of_links_f64);
        self.nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force)/inverse_temperature_in_inverse_joules
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
    pub link_length_in_meters: f64
}
impl Legendre
{
    pub fn init(number_of_links: u16) -> Legendre
    {
        Legendre
        {
            number_of_links: number_of_links,
            number_of_links_f64: number_of_links as f64,
            link_length_in_meters: 1.0
        }
    }
    // pub fn helmholtz_free_energy_in_joules<T>(&self, nondimensional_force: &T)
    // {

    // }
    // pub fn helmholtz_free_energy_per_link_in_joules<T>(&self, nondimensional_force: &T)
    // {

    // }
    pub fn relative_helmholtz_free_energy_in_joules<T>(&self, force_in_newtons: &T, inverse_temperature_in_inverse_joules: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.relative_helmholtz_free_energy_per_link_in_joules(force_in_newtons, inverse_temperature_in_inverse_joules)*self.number_of_links_f64
    }
    pub fn relative_helmholtz_free_energy_per_link_in_joules<T>(&self, force_in_newtons: &T, inverse_temperature_in_inverse_joules: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        let nondimensional_force = &(*force_in_newtons*inverse_temperature_in_inverse_joules*self.number_of_links_f64);
        self.nondimensional_relative_helmholtz_free_energy(nondimensional_force)/inverse_temperature_in_inverse_joules
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