mod test;
use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::ZERO;
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    number_of_links_f64: f64,
    contour_length: f64
}
impl FJC
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
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
    pub fn force(&self, potential_distance: &f64, potential_stiffness: &f64) -> f64
    {
        potential_stiffness*potential_distance
    }
    pub fn nondimensional_force(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64
    }
    pub fn end_to_end_length(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = potential_stiffness*potential_distance*self.link_length/BOLTZMANN_CONSTANT/temperature;
        self.contour_length*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force) - potential_stiffness*self.contour_length.powi(2)*self.link_length/BOLTZMANN_CONSTANT/temperature*((1.0/nondimensional_force.tanh() - 1.0/nondimensional_force)*(nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2)) + ((nondimensional_force.sinh()).powi(-2)/nondimensional_force.tanh() - nondimensional_force.powi(-3))/self.number_of_links_f64)
    }
    pub fn end_to_end_length_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = potential_stiffness*potential_distance*self.link_length/BOLTZMANN_CONSTANT/temperature;
        self.link_length*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force) - potential_stiffness*self.number_of_links_f64*self.link_length.powi(3)/BOLTZMANN_CONSTANT/temperature*((1.0/nondimensional_force.tanh() - 1.0/nondimensional_force)*(nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2)) + ((nondimensional_force.sinh()).powi(-2)/nondimensional_force.tanh() - nondimensional_force.powi(-3))/self.number_of_links_f64)
    }
    pub fn nondimensional_end_to_end_length(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        let nondimensional_force = nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64;
        self.number_of_links_f64*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force) - nondimensional_potential_stiffness*((1.0/nondimensional_force.tanh() - 1.0/nondimensional_force)*(nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2)) + ((nondimensional_force.sinh()).powi(-2)/nondimensional_force.tanh() - nondimensional_force.powi(-3))/self.number_of_links_f64)
    }
    pub fn nondimensional_end_to_end_length_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        let nondimensional_force = nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64;
        1.0/nondimensional_force.tanh() - 1.0/nondimensional_force - nondimensional_potential_stiffness*((1.0/nondimensional_force.tanh() - 1.0/nondimensional_force)*(nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2)) + ((nondimensional_force.sinh()).powi(-2)/nondimensional_force.tanh() - nondimensional_force.powi(-3))/self.number_of_links_f64)/self.number_of_links_f64
    }
    pub fn gibbs_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = potential_stiffness*potential_distance*self.link_length/BOLTZMANN_CONSTANT/temperature;
        -self.number_of_links_f64*BOLTZMANN_CONSTANT*temperature*(nondimensional_force.sinh()/nondimensional_force).ln() + 0.5*potential_stiffness*self.contour_length.powi(2)*((1.0/nondimensional_force.tanh() - 1.0/nondimensional_force).powi(2) + (1.0 + nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2))/self.number_of_links_f64) - self.number_of_links_f64*BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn gibbs_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = potential_stiffness*potential_distance*self.link_length/BOLTZMANN_CONSTANT/temperature;
        -BOLTZMANN_CONSTANT*temperature*(nondimensional_force.sinh()/nondimensional_force).ln() + 0.5*potential_stiffness*self.contour_length.powi(2)*((1.0/nondimensional_force.tanh() - 1.0/nondimensional_force).powi(2) + (1.0 + nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2))/self.number_of_links_f64)/self.number_of_links_f64 - BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn relative_gibbs_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.gibbs_free_energy(potential_distance, potential_stiffness, temperature) - self.gibbs_free_energy(&(ZERO*self.number_of_links_f64*self.link_length), potential_stiffness, temperature)
    }
    pub fn relative_gibbs_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.gibbs_free_energy_per_link(potential_distance, potential_stiffness, temperature) - self.gibbs_free_energy_per_link(&(ZERO*self.number_of_links_f64*self.link_length), potential_stiffness, temperature)
    }
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64;
        -self.number_of_links_f64*(nondimensional_force.sinh()/nondimensional_force).ln() + 0.5*nondimensional_potential_stiffness*((1.0/nondimensional_force.tanh() - 1.0/nondimensional_force).powi(2) + (1.0 + nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2))/self.number_of_links_f64) - self.number_of_links_f64*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64;
        -(nondimensional_force.sinh()/nondimensional_force).ln() + 0.5*nondimensional_potential_stiffness*((1.0/nondimensional_force.tanh() - 1.0/nondimensional_force).powi(2) + (1.0 + nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2))/self.number_of_links_f64)/self.number_of_links_f64 - (8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        self.nondimensional_gibbs_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness, &300.0) - self.nondimensional_gibbs_free_energy(&ZERO, nondimensional_potential_stiffness, &300.0)
    }
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        self.nondimensional_gibbs_free_energy_per_link(nondimensional_potential_distance, nondimensional_potential_stiffness, &300.0) - self.nondimensional_gibbs_free_energy_per_link(&ZERO, nondimensional_potential_stiffness, &300.0)
    }
}
