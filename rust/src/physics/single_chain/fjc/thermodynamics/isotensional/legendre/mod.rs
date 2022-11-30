pub mod test;
use std::f64::consts::PI;
use crate::math::
{
    ln,
    ln_sinhc,
    langevin
};
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::fjc::thermodynamics::IsotensionalLegendre;
pub static ZERO: f64 = 1e-6;
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub number_of_links_f64: f64,
    pub contour_length: f64
}
impl IsotensionalLegendre for FJC
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> FJC
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length
        }
    }
    fn helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        self.helmholtz_free_energy_per_link(force, temperature)*self.number_of_links_f64
    }
    fn helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy_per_link(&(*force/BOLTZMANN_CONSTANT/temperature*self.link_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    fn relative_helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        self.relative_helmholtz_free_energy_per_link(force, temperature)*self.number_of_links_f64
    }
    fn relative_helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(&(*force/BOLTZMANN_CONSTANT/temperature*self.link_length))*BOLTZMANN_CONSTANT*temperature
    }
    fn nondimensional_helmholtz_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy_per_link(nondimensional_force, temperature)*self.number_of_links_f64
    }
    fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_force) - ln(&(*nondimensional_force*0.0 + 8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)))
    }
    fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_force: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_force)*self.number_of_links_f64
    }
    fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64) -> f64
    {
        *nondimensional_force*langevin(nondimensional_force) - ln_sinhc(nondimensional_force)
    }
}