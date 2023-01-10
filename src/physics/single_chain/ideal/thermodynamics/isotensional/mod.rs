mod test;
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
    number_of_links_f64: f64
}
impl Ideal
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        Ideal
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_links_f64: number_of_links as f64
        }
    }
    pub fn end_to_end_length(&self, force: &f64, temperature: &f64) -> f64
    {
        self.number_of_links_f64*force/BOLTZMANN_CONSTANT/temperature*self.link_length.powi(2)/3.0
    }
    pub fn end_to_end_length_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        force/BOLTZMANN_CONSTANT/temperature*self.link_length.powi(2)/3.0
    }
    pub fn nondimensional_end_to_end_length(&self, nondimensional_force: &f64) -> f64
    {
        self.number_of_links_f64*nondimensional_force/3.0
    }
    pub fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_force/3.0
    }
    pub fn gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        -self.number_of_links_f64*((force*self.link_length).powi(2)/6.0/BOLTZMANN_CONSTANT/temperature + BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())
    }
    pub fn gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        -(force*self.link_length).powi(2)/6.0/BOLTZMANN_CONSTANT/temperature - BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn relative_gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        -self.number_of_links_f64*(force*self.link_length).powi(2)/6.0/BOLTZMANN_CONSTANT/temperature
    }
    pub fn relative_gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        -(force*self.link_length).powi(2)/6.0/BOLTZMANN_CONSTANT/temperature
    }
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        -self.number_of_links_f64*nondimensional_force.powi(2)/6.0 - self.number_of_links_f64*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        -nondimensional_force.powi(2)/6.0 - (8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: &f64) -> f64
    {
        -self.number_of_links_f64*nondimensional_force.powi(2)/6.0
    }
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: &f64) -> f64
    {
        -nondimensional_force.powi(2)/6.0
    }
}
