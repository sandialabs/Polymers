#[cfg(feature = "extern")]
pub mod ex;

#[cfg(feature = "python")]
pub mod py;

mod test;

use crate::physics::
{
    BOLTZMANN_CONSTANT,
    single_chain::ZERO
};
use super::
{
    nondimensional_end_to_end_length,
    nondimensional_gibbs_free_energy
};

/// The structure of the thermodynamics of the WLC model in the isotensional ensemble approximated using a Legendre transformation.
pub struct WLC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The persistance length of the chain in units of nm.
    pub persistance_length: f64,

    nondimensional_persistance_length: f64
}

/// The Helmholtz free energy as a function of the applied force and temperature, parameterized by the number of links, link length, hinge mass, and persistance length.
pub fn helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, persistance_length: &f64, force: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, &(persistance_length/((*number_of_links as f64)*link_length)), &(force*link_length/BOLTZMANN_CONSTANT/temperature), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The Helmholtz free energy per link as a function of the applied force and temperature, parameterized by the number of links, link length, hinge mass, and persistance length.
pub fn helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, persistance_length: &f64, force: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy_per_link(number_of_links, link_length, hinge_mass, &(persistance_length/((*number_of_links as f64)*link_length)), &(force*link_length/BOLTZMANN_CONSTANT/temperature), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The relative Helmholtz free energy as a function of the applied force and temperature, parameterized by the number of links, link length, and persistance length.
pub fn relative_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, persistance_length: &f64, force: &f64, temperature: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy(number_of_links, &(persistance_length/((*number_of_links as f64)*link_length)), &(force*link_length/BOLTZMANN_CONSTANT/temperature))*BOLTZMANN_CONSTANT*temperature
}

/// The relative Helmholtz free energy per link as a function of the applied force and temperature, parameterized by the number of links, link length, and persistance length.
pub fn relative_helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, persistance_length: &f64, force: &f64, temperature: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy_per_link(number_of_links, &(persistance_length/((*number_of_links as f64)*link_length)), &(force*link_length/BOLTZMANN_CONSTANT/temperature))*BOLTZMANN_CONSTANT*temperature
}

/// The nondimensional Helmholtz free energy as a function of the applied nondimensional force and temperature, parameterized by the number of links, link length, hinge mass, and nondimensional persistance length.
pub fn nondimensional_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_persistance_length: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    nondimensional_gibbs_free_energy(number_of_links, link_length, hinge_mass, nondimensional_persistance_length, nondimensional_force, temperature) + nondimensional_force*nondimensional_end_to_end_length(number_of_links, nondimensional_persistance_length, nondimensional_force)
}

/// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional force and temperature, parameterized by the number of links, link length hinge mass, and nondimensional persistance length.
pub fn nondimensional_helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_persistance_length: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, nondimensional_persistance_length, nondimensional_force, temperature)/(*number_of_links as f64)
}

/// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional force, parameterized by the number of links and nondimensional persistance length.
pub fn nondimensional_relative_helmholtz_free_energy(number_of_links: &u8, nondimensional_persistance_length: &f64, nondimensional_force: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, &1.0, &1.0, nondimensional_persistance_length, nondimensional_force, &300.0) - nondimensional_helmholtz_free_energy(number_of_links, &1.0, &1.0, nondimensional_persistance_length, &ZERO, &300.0)
}

/// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional force, parameterized by the number of links and nondimensional persistance length.
pub fn nondimensional_relative_helmholtz_free_energy_per_link(number_of_links: &u8, nondimensional_persistance_length: &f64, nondimensional_force: &f64) -> f64
{
    nondimensional_helmholtz_free_energy_per_link(number_of_links, &1.0, &1.0, nondimensional_persistance_length, nondimensional_force, &300.0) - nondimensional_helmholtz_free_energy_per_link(number_of_links, &1.0, &1.0, nondimensional_persistance_length, &ZERO, &300.0)
}

/// The implemented functionality of the thermodynamics of the WLC model in the isotensional ensemble approximated using a Legendre transformation.
impl WLC
{
    /// Initializes and returns an instance of the thermodynamics of the WLC model in the isotensional ensemble approximated using a Legendre transformation.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, persistance_length: f64) -> Self
    {
        WLC
        {
            hinge_mass,
            link_length,
            number_of_links,
            persistance_length,
            nondimensional_persistance_length: persistance_length/(number_of_links as f64)/link_length
        }
    }
    /// The Helmholtz free energy as a function of the applied force and temperature.
    pub fn helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.persistance_length, force, temperature)
    }
    /// The Helmholtz free energy per link as a function of the applied force and temperature.
    pub fn helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.persistance_length, force, temperature)
    }
    /// The relative Helmholtz free energy as a function of the applied force and temperature.
    pub fn relative_helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.persistance_length, force, temperature)
    }
    /// The relative Helmholtz free energy per link as a function of the applied force and temperature.
    pub fn relative_helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.persistance_length, force, temperature)
    }
    /// The nondimensional Helmholtz free energy as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.nondimensional_persistance_length, nondimensional_force, temperature)
    }
    /// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.nondimensional_persistance_length, nondimensional_force, temperature)
    }
    /// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional force.
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy(&self.number_of_links, &self.nondimensional_persistance_length, nondimensional_force)
    }
    /// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional force.
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy_per_link(&self.number_of_links, &self.nondimensional_persistance_length, nondimensional_force)
    }
}
