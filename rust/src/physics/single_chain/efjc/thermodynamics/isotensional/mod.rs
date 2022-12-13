pub mod test;
pub mod asymptotic;
pub mod legendre;
use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::fjc::ZERO;
pub struct EFJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub link_stiffness: f64,
    pub number_of_links_f64: f64,
    pub contour_length: f64,
    pub legendre: self::legendre::EFJC
}
use super::Isotensional;
impl Isotensional for EFJC
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> EFJC
    {
        EFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length,
            legendre: self::legendre::EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness)
        }
    }
}
pub trait Legendre
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self;
//    fn helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64;
//    fn helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64;
//    fn relative_helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64;
//    fn relative_helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64;
//    fn nondimensional_helmholtz_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
//    fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
//    fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_force: &f64) -> f64;
//    fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64) -> f64;
}
pub trait Asymptotic
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self;
    fn end_to_end_length(&self, force: &f64, temperature: &f64) -> f64;
    fn end_to_end_length_per_link(&self, force: &f64, temperature: &f64) -> f64;
    fn nondimensional_end_to_end_length(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
    fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
    fn gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64;
    fn gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64;
    fn relative_gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64;
    fn relative_gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64;
    fn nondimensional_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
    fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
    fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
    fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
}
