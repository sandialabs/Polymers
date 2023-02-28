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

/// The structure of the thermodynamics of the FJC model in the isotensional ensemble approximated using a Legendre transformation.
pub struct FJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8
}

/// The Helmholtz free energy as a function of the applied force and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, force: &f64, temperature: &f64) -> f64
{
    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
    (*number_of_links as f64)*BOLTZMANN_CONSTANT*temperature*(nondimensional_force/nondimensional_force.tanh() - 1.0 - (nondimensional_force.sinh()/nondimensional_force).ln() - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())
}

/// The Helmholtz free energy per link as a function of the applied force and temperature, parameterized by the link length and hinge mass.
pub fn helmholtz_free_energy_per_link(link_length: &f64, hinge_mass: &f64, force: &f64, temperature: &f64) -> f64
{
    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
    BOLTZMANN_CONSTANT*temperature*(nondimensional_force/nondimensional_force.tanh() - 1.0 - (nondimensional_force.sinh()/nondimensional_force).ln() - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())
}

/// The relative Helmholtz free energy as a function of the applied force and temperature, parameterized by the number of links and link length.
pub fn relative_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, force: &f64, temperature: &f64) -> f64
{
    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
    (*number_of_links as f64)*BOLTZMANN_CONSTANT*temperature*(nondimensional_force/nondimensional_force.tanh() - 1.0 - (nondimensional_force.sinh()/nondimensional_force).ln())
}

/// The relative Helmholtz free energy per link as a function of the applied force and temperature, parameterized by the link length.
pub fn relative_helmholtz_free_energy_per_link(link_length: &f64, force: &f64, temperature: &f64) -> f64
{
    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
    BOLTZMANN_CONSTANT*temperature*(nondimensional_force/nondimensional_force.tanh() - 1.0 - (nondimensional_force.sinh()/nondimensional_force).ln())
}

/// The nondimensional Helmholtz free energy as a function of the applied nondimensional force and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    (*number_of_links as f64)*(nondimensional_force/nondimensional_force.tanh() - 1.0 - (nondimensional_force.sinh()/nondimensional_force).ln() - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())
}

/// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional force and temperature, parameterized by the link length and hinge mass.
pub fn nondimensional_helmholtz_free_energy_per_link(link_length: &f64, hinge_mass: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    nondimensional_force/nondimensional_force.tanh() - 1.0 - (nondimensional_force.sinh()/nondimensional_force).ln() - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional force, parameterized by the number of links.
pub fn nondimensional_relative_helmholtz_free_energy(number_of_links: &u8, nondimensional_force: &f64) -> f64
{
    (*number_of_links as f64)*(nondimensional_force/nondimensional_force.tanh() - 1.0 - (nondimensional_force.sinh()/nondimensional_force).ln())
}

/// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional force.
pub fn nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_force: &f64) -> f64
{
    nondimensional_force/nondimensional_force.tanh() - 1.0 - (nondimensional_force.sinh()/nondimensional_force).ln()
}

/// The implemented functionality of the thermodynamics of the FJC model in the isotensional ensemble approximated using a Legendre transformation.
impl FJC
{
    /// Initializes and returns an instance of the thermodynamics of the FJC model in the isotensional ensemble approximated using a Legendre transformation.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links
        }
    }
    /// The Helmholtz free energy as a function of the applied force and temperature.
    pub fn helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, force, temperature)
    }
    /// The Helmholtz free energy per link as a function of the applied force and temperature.
    pub fn helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy_per_link(&self.link_length, &self.hinge_mass, force, temperature)
    }
    /// The relative Helmholtz free energy as a function of the applied force and temperature.
    pub fn relative_helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy(&self.number_of_links, &self.link_length, force, temperature)
    }
    /// The relative Helmholtz free energy per link as a function of the applied force and temperature.
    pub fn relative_helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy_per_link(&self.link_length, force, temperature)
    }
    /// The nondimensional Helmholtz free energy as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, nondimensional_force, temperature)
    }
    /// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy_per_link(&self.link_length, &self.hinge_mass, nondimensional_force, temperature)
    }
    /// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional force.
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy(&self.number_of_links, nondimensional_force)
    }
    /// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional force.
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_force)
    }
}
