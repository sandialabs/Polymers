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

/// The structure of the thermodynamics of the ideal chain model in the isotensional ensemble.
pub struct Ideal
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8
}

/// The expected end-to-end length as a function of the applied force and temperature, parameterized by the number of links and link length.
pub fn end_to_end_length(number_of_links: &u8, link_length: &f64, force: &f64, temperature: &f64) -> f64
{
    (*number_of_links as f64)*force/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/3.0
}

/// The expected end-to-end length per link as a function of the applied force and temperature, parameterized by the number of links and link length.
pub fn end_to_end_length_per_link(link_length: &f64, force: &f64, temperature: &f64) -> f64
{
    force/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/3.0
}

/// The expected nondimensional end-to-end length as a function of the applied nondimensional force, parameterized by the number of links.
pub fn nondimensional_end_to_end_length(number_of_links: &u8, nondimensional_force: &f64) -> f64
{
    (*number_of_links as f64)*nondimensional_force/3.0
}

/// The expected nondimensional end-to-end length per link as a function of the applied nondimensional force.
pub fn nondimensional_end_to_end_length_per_link(nondimensional_force: &f64) -> f64
{
    nondimensional_force/3.0
}

/// The Gibbs free energy as a function of the applied nondimensional force and temperature, parameterized by the number of links and link length.
pub fn gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, force: &f64, temperature: &f64) -> f64
{
    -(*number_of_links as f64)*((force*link_length).powi(2)/6.0/BOLTZMANN_CONSTANT/temperature + BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())
}

/// The Gibbs free energy per link as a function of the applied nondimensional force and temperature, parameterized by the link length.
pub fn gibbs_free_energy_per_link(link_length: &f64, hinge_mass: &f64, force: &f64, temperature: &f64) -> f64
{
    -(force*link_length).powi(2)/6.0/BOLTZMANN_CONSTANT/temperature - BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The relative Gibbs free energy as a function of the applied nondimensional force and temperature, parameterized by the number of links and link length.
pub fn relative_gibbs_free_energy(number_of_links: &u8, link_length: &f64, force: &f64, temperature: &f64) -> f64
{
    -(*number_of_links as f64)*(force*link_length).powi(2)/6.0/BOLTZMANN_CONSTANT/temperature
}

/// The relative Gibbs free energy per link as a function of the applied nondimensional force and temperature, parameterized by the link length.
pub fn relative_gibbs_free_energy_per_link(link_length: &f64, force: &f64, temperature: &f64) -> f64
{
    -(force*link_length).powi(2)/6.0/BOLTZMANN_CONSTANT/temperature
}

/// The nondimensional Gibbs free energy as a function of the applied nondimensional force and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    -(*number_of_links as f64)*(nondimensional_force.powi(2)/6.0 + (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())
}

/// The nondimensional Gibbs free energy per link as a function of the applied nondimensional force and temperature, parameterized by the link length and hinge mass.
pub fn nondimensional_gibbs_free_energy_per_link(link_length: &f64, hinge_mass: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    -nondimensional_force.powi(2)/6.0 - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional relative Gibbs free energy as a function of the applied nondimensional force, parameterized by the number of links.
pub fn nondimensional_relative_gibbs_free_energy(number_of_links: &u8, nondimensional_force: &f64) -> f64
{
    -(*number_of_links as f64)*nondimensional_force.powi(2)/6.0
}

/// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional force.
pub fn nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force: &f64) -> f64
{
    -nondimensional_force.powi(2)/6.0
}

/// The implemented functionality of the thermodynamics of the ideal chain model in the isotensional ensemble.
impl Ideal
{
    /// Initializes and returns an instance of the thermodynamics of the ideal chain model in the isotensional ensemble.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        Ideal
        {
            hinge_mass,
            link_length,
            number_of_links
        }
    }
    /// The expected end-to-end length as a function of the applied force and temperature.
    pub fn end_to_end_length(&self, force: &f64, temperature: &f64) -> f64
    {
        end_to_end_length(&self.number_of_links, &self.link_length, force, temperature)
    }
    /// The expected end-to-end length per link as a function of the applied force and temperature.
    pub fn end_to_end_length_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        end_to_end_length_per_link(&self.link_length, force, temperature)
    }
    /// The expected nondimensional end-to-end length as a function of the applied nondimensional force.
    pub fn nondimensional_end_to_end_length(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_end_to_end_length(&self.number_of_links, nondimensional_force)
    }
    /// The expected nondimensional end-to-end length per link as a function of the applied nondimensional force.
    pub fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_end_to_end_length_per_link(nondimensional_force)
    }
    /// The Gibbs free energy as a function of the applied force and temperature.
    pub fn gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, force, temperature)
    }
    /// The Gibbs free energy per link as a function of the applied force and temperature.
    pub fn gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        gibbs_free_energy_per_link(&self.link_length, &self.hinge_mass, force, temperature)
    }
    /// The relative Gibbs free energy as a function of the applied force and temperature.
    pub fn relative_gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_gibbs_free_energy(&self.number_of_links, &self.link_length, force, temperature)
    }
    /// The relative Gibbs free energy per link as a function of the applied force and temperature.
    pub fn relative_gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_gibbs_free_energy_per_link(&self.link_length, force, temperature)
    }
    /// The nondimensional Gibbs free energy as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, nondimensional_force, temperature)
    }
    /// The nondimensional Gibbs free energy per link as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_gibbs_free_energy_per_link(&self.link_length, &self.hinge_mass, nondimensional_force, temperature)
    }
    /// The nondimensional relative Gibbs free energy as a function of the applied nondimensional force.
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_relative_gibbs_free_energy(&self.number_of_links, nondimensional_force)
    }
    /// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional force.
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force)
    }
}
