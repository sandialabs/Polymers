pub mod test;
use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
pub struct Ideal
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub number_of_links_f64: f64,
    pub contour_length: f64
}
use super::Isotensional;
impl Isotensional for Ideal
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Ideal
    {
        Ideal
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length
        }
    }
    fn end_to_end_length(&self, force: &f64, temperature: &f64) -> f64
    {
        self.number_of_links_f64*force/BOLTZMANN_CONSTANT/temperature*self.link_length.powi(2)/3.0
    }
    fn end_to_end_length_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        force/BOLTZMANN_CONSTANT/temperature*self.link_length.powi(2)/3.0
    }
    fn nondimensional_end_to_end_length(&self, nondimensional_force: &f64) -> f64
    {
        self.number_of_links_f64*nondimensional_force/3.0
    }
    fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_force/3.0
    }
    fn gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        -self.number_of_links_f64*((force*self.link_length).powi(2)/6.0/BOLTZMANN_CONSTANT/temperature + BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())
    }
    fn gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        -(force*self.link_length).powi(2)/6.0/BOLTZMANN_CONSTANT/temperature - BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    fn relative_gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        -self.number_of_links_f64*(force*self.link_length).powi(2)/6.0/BOLTZMANN_CONSTANT/temperature
    }
    fn relative_gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        -(force*self.link_length).powi(2)/6.0/BOLTZMANN_CONSTANT/temperature
    }
    fn nondimensional_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        -self.number_of_links_f64*nondimensional_force.powi(2)/6.0 - self.number_of_links_f64*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        -nondimensional_force.powi(2)/6.0 - (8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: &f64) -> f64
    {
        -self.number_of_links_f64*nondimensional_force.powi(2)/6.0
    }
    fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: &f64) -> f64
    {
        -nondimensional_force.powi(2)/6.0
    }
}