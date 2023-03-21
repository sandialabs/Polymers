#[cfg(feature = "extern")]
pub mod ex;

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
    pub number_of_links: u8
}

/// The expected end-to-end length as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links and link length.
pub fn end_to_end_length(number_of_links: &u8, link_length: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    let nondimensional_force = potential_stiffness*potential_distance*link_length/BOLTZMANN_CONSTANT/temperature;
    (*number_of_links as f64)*link_length*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force)*(1.0 - potential_stiffness*(*number_of_links as f64)*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature*(nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2)))
}

/// The expected end-to-end length per link as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links and link length.
pub fn end_to_end_length_per_link(number_of_links: &u8, link_length: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    let nondimensional_force = potential_stiffness*potential_distance*link_length/BOLTZMANN_CONSTANT/temperature;
    link_length*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force)*(1.0 - potential_stiffness*(*number_of_links as f64)*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature*(nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2)))
}

/// The expected nondimensional end-to-end length as a function of the applied nondimensional potential distance and nondimensional potential stiffness, parameterized by the number of links.
pub fn nondimensional_end_to_end_length(number_of_links: &u8, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
{
    let nondimensional_force = (*number_of_links as f64)*nondimensional_potential_stiffness*nondimensional_potential_distance;
    (*number_of_links as f64)*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force)*(1.0 - (*number_of_links as f64)*nondimensional_potential_stiffness*(nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2)))
}

/// The expected nondimensional end-to-end length per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness, parameterized by the number of links.
pub fn nondimensional_end_to_end_length_per_link(number_of_links: &u8, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
{
    let nondimensional_force = (*number_of_links as f64)*nondimensional_potential_stiffness*nondimensional_potential_distance;
    (1.0/nondimensional_force.tanh() - 1.0/nondimensional_force)*(1.0 - (*number_of_links as f64)*nondimensional_potential_stiffness*(nondimensional_force.powi(-2) - (nondimensional_force.sinh()).powi(-2)))
}

/// The expected force as a function of the applied potential distance, potential stiffness, and temperature.
pub fn force(potential_distance: &f64, potential_stiffness: &f64) -> f64
{
    potential_stiffness*potential_distance
}

/// The expected nondimensional force as a function of the applied nondimensional potential distance and nondimensional potential stiffness, parameterized by the number of links.
pub fn nondimensional_force(number_of_links: &u8, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
{
    (*number_of_links as f64)*nondimensional_potential_stiffness*nondimensional_potential_distance
}

/// The Gibbs free energy as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    let nondimensional_force = potential_stiffness*potential_distance*link_length/BOLTZMANN_CONSTANT/temperature;
    -(*number_of_links as f64)*BOLTZMANN_CONSTANT*temperature*(nondimensional_force.sinh()/nondimensional_force).ln() + 0.5*potential_stiffness*((*number_of_links as f64)*link_length).powi(2)*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force).powi(2) - (*number_of_links as f64)*BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The Gibbs free energy epr link as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    let nondimensional_force = potential_stiffness*potential_distance*link_length/BOLTZMANN_CONSTANT/temperature;
    -BOLTZMANN_CONSTANT*temperature*(nondimensional_force.sinh()/nondimensional_force).ln() + 0.5*potential_stiffness*((*number_of_links as f64)*link_length).powi(2)/(*number_of_links as f64)*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force).powi(2) - BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The relative Gibbs free energy as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links and link length.
pub fn relative_gibbs_free_energy(number_of_links: &u8, link_length: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    gibbs_free_energy(number_of_links, link_length, &1.0, potential_distance, potential_stiffness, temperature) - gibbs_free_energy(number_of_links, link_length, &1.0, &(ZERO*(*number_of_links as f64)*link_length), potential_stiffness, temperature)
}

/// The relative Gibbs free energy per link as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links and link length.
pub fn relative_gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    gibbs_free_energy_per_link(number_of_links, link_length, &1.0, potential_distance, potential_stiffness, temperature) - gibbs_free_energy_per_link(number_of_links, link_length, &1.0, &(ZERO*(*number_of_links as f64)*link_length), potential_stiffness, temperature)
}

/// The nondimensional Gibbs free energy as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
{
    let nondimensional_force = (*number_of_links as f64)*nondimensional_potential_stiffness*nondimensional_potential_distance;
    -(*number_of_links as f64)*(nondimensional_force.sinh()/nondimensional_force).ln() + 0.5*nondimensional_potential_stiffness*(*number_of_links as f64).powi(2)*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force).powi(2) - (*number_of_links as f64)*(8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional Gibbs free energy per link as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
{
    let nondimensional_force = (*number_of_links as f64)*nondimensional_potential_stiffness*nondimensional_potential_distance;
    -(nondimensional_force.sinh()/nondimensional_force).ln() + 0.5*nondimensional_potential_stiffness*(*number_of_links as f64)*(1.0/nondimensional_force.tanh() - 1.0/nondimensional_force).powi(2) - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional relative Gibbs free energy as a function of the applied nondimensional potential distance and nondimensional potential stiffness, parameterized by the number of links.
pub fn nondimensional_relative_gibbs_free_energy(number_of_links: &u8, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
{
    nondimensional_gibbs_free_energy(number_of_links, &1.0, &1.0, nondimensional_potential_distance, nondimensional_potential_stiffness, &300.0) - nondimensional_gibbs_free_energy(number_of_links, &1.0, &1.0, &ZERO, nondimensional_potential_stiffness, &300.0)
}

/// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness, parameterized by the number of links.
pub fn nondimensional_relative_gibbs_free_energy_per_link(number_of_links: &u8, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
{
    nondimensional_gibbs_free_energy_per_link(number_of_links, &1.0, &1.0, nondimensional_potential_distance, nondimensional_potential_stiffness, &300.0) - nondimensional_gibbs_free_energy_per_link(number_of_links, &1.0, &1.0, &ZERO, nondimensional_potential_stiffness, &300.0)
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
            number_of_links
        }
    }
    /// The expected end-to-end length as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn end_to_end_length(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        end_to_end_length(&self.number_of_links, &self.link_length, potential_distance, potential_stiffness, temperature)
    }
    /// The expected end-to-end length per link as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn end_to_end_length_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        end_to_end_length_per_link(&self.number_of_links, &self.link_length, potential_distance, potential_stiffness, temperature)
    }
    /// The expected nondimensional end-to-end length as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    pub fn nondimensional_end_to_end_length(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        nondimensional_end_to_end_length(&self.number_of_links, nondimensional_potential_distance, nondimensional_potential_stiffness)
    }
    /// The expected nondimensional end-to-end length per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    pub fn nondimensional_end_to_end_length_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        nondimensional_end_to_end_length_per_link(&self.number_of_links, nondimensional_potential_distance, nondimensional_potential_stiffness)
    }
    /// The expected force as a function of the applied potential distance and potential stiffness
    pub fn force(&self, potential_distance: &f64, potential_stiffness: &f64) -> f64
    {
        force(potential_distance, potential_stiffness)
    }
    /// The expected nondimensional force as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    pub fn nondimensional_force(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        nondimensional_force(&self.number_of_links, nondimensional_potential_distance, nondimensional_potential_stiffness)
    }
    /// The Gibbs free energy as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn gibbs_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, potential_distance, potential_stiffness, temperature)
    }
    /// The Gibbs free energy epr link as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn gibbs_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, potential_distance, potential_stiffness, temperature)
    }
    /// The relative Gibbs free energy as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn relative_gibbs_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        relative_gibbs_free_energy(&self.number_of_links, &self.link_length, potential_distance, potential_stiffness, temperature)
    }
    /// The relative Gibbs free energy per link as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn relative_gibbs_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        relative_gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, potential_distance, potential_stiffness, temperature)
    }
    /// The nondimensional Gibbs free energy as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        nondimensional_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, nondimensional_potential_distance, nondimensional_potential_stiffness, temperature)
    }
    /// The nondimensional Gibbs free energy per link as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        nondimensional_gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, nondimensional_potential_distance, nondimensional_potential_stiffness, temperature)
    }
    /// The nondimensional relative Gibbs free energy as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        nondimensional_relative_gibbs_free_energy(&self.number_of_links, nondimensional_potential_distance, nondimensional_potential_stiffness)
    }
    /// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        nondimensional_relative_gibbs_free_energy_per_link(&self.number_of_links, nondimensional_potential_distance, nondimensional_potential_stiffness)
    }
}
