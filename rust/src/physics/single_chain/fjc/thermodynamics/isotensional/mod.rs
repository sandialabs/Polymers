pub mod test;
pub mod legendre;
use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::fjc::thermodynamics::
{
    Isotensional,
    IsotensionalLegendre
};
use crate::physics::single_chain::fjc::thermodynamics::isotensional::legendre::FJC as FJCLegendre;
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub number_of_links_f64: f64,
    pub contour_length: f64,
    pub legendre: FJCLegendre
}
impl Isotensional for FJC
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> FJC
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
    fn end_to_end_length(&self, force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_end_to_end_length(&(force/BOLTZMANN_CONSTANT/temperature*self.link_length))*self.link_length
    }
    fn end_to_end_length_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_end_to_end_length_per_link(&(force/BOLTZMANN_CONSTANT/temperature*self.link_length))*self.link_length
    }
    fn nondimensional_end_to_end_length(&self, nondimensional_force: &f64) -> f64
    {
        (1.0/nondimensional_force.tanh() - 1.0/nondimensional_force)*self.number_of_links_f64
    }
    fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: &f64) -> f64
    {
        1.0/nondimensional_force.tanh() - 1.0/nondimensional_force
    }
    fn gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        self.relative_gibbs_free_energy(force, temperature) - self.number_of_links_f64*BOLTZMANN_CONSTANT*temperature*(8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)).ln()
    }
    fn gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        self.relative_gibbs_free_energy_per_link(force, temperature) - BOLTZMANN_CONSTANT*temperature*(8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)).ln()
    }
    fn relative_gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_gibbs_free_energy(&(force/BOLTZMANN_CONSTANT/temperature*self.link_length))*BOLTZMANN_CONSTANT*temperature
    }
    fn relative_gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_gibbs_free_energy_per_link(&(force/BOLTZMANN_CONSTANT/temperature*self.link_length))*BOLTZMANN_CONSTANT*temperature
    }
    fn nondimensional_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        -self.number_of_links_f64*(nondimensional_force.sinh()/nondimensional_force).ln() - self.number_of_links_f64*(8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)).ln()
    }
    fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        -(nondimensional_force.sinh()/nondimensional_force).ln() - (8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)).ln()
    }
    fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: &f64) -> f64
    {
        -self.number_of_links_f64*(nondimensional_force.sinh()/nondimensional_force).ln()
    }
    fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: &f64) -> f64
    {
        -(nondimensional_force.sinh()/nondimensional_force).ln()
    }
}
