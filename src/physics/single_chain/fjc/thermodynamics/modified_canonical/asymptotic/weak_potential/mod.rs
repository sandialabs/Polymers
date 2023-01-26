#[cfg(feature = "python")]
pub mod py;

mod test;

use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::ZERO;

/// The structure of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach valid for weak potentials.
pub struct FJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    number_of_links_f64: f64,

    contour_length: f64
}

/// The implemented functionality of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach valid for weak potentials.
impl FJC
{
    /// Initializes and returns an instance of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach valid for weak potentials.
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
    /// The expected end-to-end length as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn end_to_end_length(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = potential_stiffness*potential_distance*self.link_length/BOLTZMANN_CONSTANT/temperature;
        self.contour_length*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force)*(1.0 - potential_stiffness*self.number_of_links_f64*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature*(nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2)))
    }
    /// The expected end-to-end length per link as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn end_to_end_length_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = potential_stiffness*potential_distance*self.link_length/BOLTZMANN_CONSTANT/temperature;
        self.link_length*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force)*(1.0 - potential_stiffness*self.number_of_links_f64*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature*(nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2)))
    }
    /// The expected nondimensional end-to-end length as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    pub fn nondimensional_end_to_end_length(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        let nondimensional_force = nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64;
        self.number_of_links_f64*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force)*(1.0 - nondimensional_potential_stiffness/self.number_of_links_f64*(nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2)))
    }
    /// The expected nondimensional end-to-end length per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    pub fn nondimensional_end_to_end_length_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        let nondimensional_force = nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64;
        (1.0/nondimensional_force.tanh() - 1.0/nondimensional_force)*(1.0 - nondimensional_potential_stiffness/self.number_of_links_f64*(nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2)))
    }
    /// The expected force as a function of the applied potential distance and potential stiffness
    pub fn force(&self, potential_distance: &f64, potential_stiffness: &f64) -> f64
    {
        potential_stiffness*potential_distance
    }
    /// The expected nondimensional force as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    pub fn nondimensional_force(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64
    }
    /// The gibbs free energy as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn gibbs_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = potential_stiffness*potential_distance*self.link_length/BOLTZMANN_CONSTANT/temperature;
        -self.number_of_links_f64*BOLTZMANN_CONSTANT*temperature*(nondimensional_force.sinh()/nondimensional_force).ln() + 0.5*potential_stiffness*self.contour_length.powi(2)*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force).powi(2) - self.number_of_links_f64*BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    /// The gibbs free energy epr link as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn gibbs_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = potential_stiffness*potential_distance*self.link_length/BOLTZMANN_CONSTANT/temperature;
        -BOLTZMANN_CONSTANT*temperature*(nondimensional_force.sinh()/nondimensional_force).ln() + 0.5*potential_stiffness*self.contour_length.powi(2)/self.number_of_links_f64*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force).powi(2) - BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    /// The relative gibbs free energy as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn relative_gibbs_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.gibbs_free_energy(potential_distance, potential_stiffness, temperature) - self.gibbs_free_energy(&(ZERO*self.number_of_links_f64*self.link_length), potential_stiffness, temperature)
    }
    /// The relative gibbs free energy per link as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn relative_gibbs_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.gibbs_free_energy_per_link(potential_distance, potential_stiffness, temperature) - self.gibbs_free_energy_per_link(&(ZERO*self.number_of_links_f64*self.link_length), potential_stiffness, temperature)
    }
    /// The nondimensional gibbs free energy as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64;
        -self.number_of_links_f64*(nondimensional_force.sinh()/nondimensional_force).ln() + 0.5*nondimensional_potential_stiffness*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force).powi(2) - self.number_of_links_f64*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    /// The nondimensional gibbs free energy per link as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64;
        -(nondimensional_force.sinh()/nondimensional_force).ln() + 0.5*nondimensional_potential_stiffness/self.number_of_links_f64*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force).powi(2) - (8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    /// The nondimensional relative gibbs free energy as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        self.nondimensional_gibbs_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness, &300.0) - self.nondimensional_gibbs_free_energy(&ZERO, nondimensional_potential_stiffness, &300.0)
    }
    /// The nondimensional relative gibbs free energy per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        self.nondimensional_gibbs_free_energy_per_link(nondimensional_potential_distance, nondimensional_potential_stiffness, &300.0) - self.nondimensional_gibbs_free_energy_per_link(&ZERO, nondimensional_potential_stiffness, &300.0)
    }
}
