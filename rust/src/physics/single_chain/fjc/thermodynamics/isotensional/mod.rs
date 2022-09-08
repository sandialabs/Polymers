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
    pub legendre: Legendre
}

impl Isotensional
{
    pub fn init(number_of_links: u16) -> Isotensional
    {
        Isotensional
        {
            number_of_links: number_of_links,
            legendre: Legendre::init(number_of_links),
        }
    }
    pub fn nondimensional_end_to_end_length_per_link<T>(&self, nondimensional_force: T) -> T
    where T:
        Math<T>
    {
        langevin(&nondimensional_force)
    }
    pub fn nondimensional_relative_gibbs_free_energy_per_link<T>(&self, nondimensional_force: T) -> T
    where T:
        Math<T>,
        T: std::ops::Neg<Output = T>
    {
        -ln_sinhc(&nondimensional_force)
    }
}

pub struct Legendre
{
    pub number_of_links: u16
}

impl Legendre
{
    pub fn init(number_of_links: u16) -> Legendre
    {
        Legendre
        {
            number_of_links: number_of_links
        }
    }
    pub fn nondimensional_relative_helmholtz_free_energy_per_link<T>(&self, nondimensional_force: T) -> T
    where T:
        Math<T>,
        T: std::marker::Copy,
        T: std::ops::Mul<T, Output = T>,
        T: std::ops::Sub<T, Output = T>,
    {
        nondimensional_force*langevin(&nondimensional_force) - ln_sinhc(&nondimensional_force)
    }
}