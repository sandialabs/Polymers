use crate::math::
{
    Math,
    langevin,
    ln_sinhc
};

pub mod test;

pub struct Isometric
{
    pub number_of_links: u16
}

impl Isometric
{
    pub fn init(number_of_links: u16) -> Isometric
    {
        Isometric
        {
            number_of_links: number_of_links
        }
    }
    pub fn nondimensional_relative_helmholtz_free_energy_per_link_legendre_transformation<T>(&self, nondimensional_force: T) -> T
    where T:
        Math<T>,
        T: std::marker::Copy,
        T: std::ops::Mul<T, Output = T>,
        T: std::ops::Sub<T, Output = T>,
    {
        langevin(&nondimensional_force)*nondimensional_force - ln_sinhc(&nondimensional_force)
    }
}