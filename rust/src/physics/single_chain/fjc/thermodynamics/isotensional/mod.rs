use std::f64::consts::PI;
use crate::math::
{
    Math,
    ln,
    ln_sinhc,
    langevin
};
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::
{
    Isotensional,
    IsotensionalLegendre
};

pub mod test;

pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u16,
    pub number_of_links_f64: f64,
    pub contour_length: f64,
    pub legendre: FJCLegendre
}

impl Isotensional for FJC
{
    fn init(number_of_links: u16, link_length: f64, hinge_mass: f64) -> FJC
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length,
            legendre: FJCLegendre::init(number_of_links, link_length, hinge_mass)
        }
    }
    fn end_to_end_length<T>(&self, force: &T, temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.end_to_end_length_per_link(force, temperature)*self.number_of_links_f64
    }
    fn end_to_end_length_per_link<T>(&self, force: &T, temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_end_to_end_length_per_link(&(*force/BOLTZMANN_CONSTANT/temperature*self.link_length))*self.link_length
    }
    fn nondimensional_end_to_end_length<T>(&self, nondimensional_force: &T) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
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
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        langevin(nondimensional_force)
    }
    fn gibbs_free_energy<T>(&self, force: &T, temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.gibbs_free_energy_per_link(force, temperature)*self.number_of_links_f64
    }
    fn gibbs_free_energy_per_link<T>(&self, force: &T, temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_gibbs_free_energy_per_link(&(*force/BOLTZMANN_CONSTANT/temperature*self.link_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    fn relative_gibbs_free_energy<T>(&self, force: &T, temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.relative_gibbs_free_energy_per_link(force, temperature)*self.number_of_links_f64
    }
    fn relative_gibbs_free_energy_per_link<T>(&self, force: &T, temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_relative_gibbs_free_energy_per_link(&(*force/BOLTZMANN_CONSTANT/temperature*self.link_length))*BOLTZMANN_CONSTANT*temperature
    }
    fn nondimensional_gibbs_free_energy<T>(&self, nondimensional_force: &T, temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_gibbs_free_energy_per_link(nondimensional_force, temperature)*self.number_of_links_f64
    }
    fn nondimensional_gibbs_free_energy_per_link<T>(&self, nondimensional_force: &T, temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force) - ln(&(*nondimensional_force*0.0 + 8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)))
    }
    fn nondimensional_relative_gibbs_free_energy<T>(&self, nondimensional_force: &T) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
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
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        -ln_sinhc(nondimensional_force)
    }
}

pub struct FJCLegendre
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u16,
    pub number_of_links_f64: f64,
    pub contour_length: f64
}

impl IsotensionalLegendre for FJCLegendre
{
    fn init(number_of_links: u16, link_length: f64, hinge_mass: f64) -> FJCLegendre
    {
        FJCLegendre
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length
        }
    }
    fn helmholtz_free_energy<T>(&self, force: &T, temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.helmholtz_free_energy_per_link(force, temperature)*self.number_of_links_f64
    }
    fn helmholtz_free_energy_per_link<T>(&self, force: &T, temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_helmholtz_free_energy_per_link(&(*force/BOLTZMANN_CONSTANT/temperature*self.link_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    fn relative_helmholtz_free_energy<T>(&self, force: &T, temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.relative_helmholtz_free_energy_per_link(force, temperature)*self.number_of_links_f64
    }
    fn relative_helmholtz_free_energy_per_link<T>(&self, force: &T, temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(&(*force/BOLTZMANN_CONSTANT/temperature*self.link_length))*BOLTZMANN_CONSTANT*temperature
    }
    fn nondimensional_helmholtz_free_energy<T>(&self, nondimensional_force: &T, temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_helmholtz_free_energy_per_link(nondimensional_force, temperature)*self.number_of_links_f64
    }
    fn nondimensional_helmholtz_free_energy_per_link<T>(&self, nondimensional_force: &T, temperature: f64) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_force) - ln(&(*nondimensional_force*0.0 + 8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)))
    }
    fn nondimensional_relative_helmholtz_free_energy<T>(&self, nondimensional_force: &T) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_force)*self.number_of_links_f64
    }
    fn nondimensional_relative_helmholtz_free_energy_per_link<T>(&self, nondimensional_force: &T) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        *nondimensional_force*langevin(nondimensional_force) - ln_sinhc(nondimensional_force)
    }
}
