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

/// The structure of the thermodynamics of the SWFJC model in the isotensional ensemble approximated using a Legendre transformation.
pub struct SWFJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The width of the well in units of nm.
    pub well_width: f64
}

/// The Helmholtz free energy as a function of the applied force and temperature, parameterized by the number of links, link length, hinge mass, and well width.
pub fn helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, well_width: &f64, force: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, well_width, &(*force/BOLTZMANN_CONSTANT/temperature*link_length), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The Helmholtz free energy per link as a function of the applied force and temperature, parameterized by the link length, hinge mass, and well width.
pub fn helmholtz_free_energy_per_link(link_length: &f64, hinge_mass: &f64, well_width: &f64, force: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy_per_link(link_length, hinge_mass, well_width, &(*force/BOLTZMANN_CONSTANT/temperature*link_length), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The relative Helmholtz free energy as a function of the applied force and temperature, parameterized by the number of links, link length, and well width.
pub fn relative_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, well_width: &f64, force: &f64, temperature: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy(number_of_links, link_length, well_width, &(*force/BOLTZMANN_CONSTANT/temperature*link_length))*BOLTZMANN_CONSTANT*temperature
}

/// The relative Helmholtz free energy per link as a function of the applied force and temperature, parameterized by the link length and well width.
pub fn relative_helmholtz_free_energy_per_link(link_length: &f64, well_width: &f64, force: &f64, temperature: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy_per_link(link_length, well_width, &(*force/BOLTZMANN_CONSTANT/temperature*link_length))*BOLTZMANN_CONSTANT*temperature
}

/// The nondimensional Helmholtz free energy as a function of the applied nondimensional force and temperature, parameterized by the number of links, link length, hinge mass, and well width.
pub fn nondimensional_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, well_width: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    (*number_of_links as f64)*nondimensional_helmholtz_free_energy_per_link(&link_length, &hinge_mass, &well_width, nondimensional_force, temperature)
}

/// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional force and temperature, parameterized by the link length, hinge mass, and well width.
pub fn nondimensional_helmholtz_free_energy_per_link(link_length: &f64, hinge_mass: &f64, well_width: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    let nondimensional_well_parameter = 1.0 + well_width/link_length;
    (nondimensional_well_parameter.powi(2)*nondimensional_force.powi(2)*(nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force.powi(2)*nondimensional_force.sinh())/(nondimensional_well_parameter*nondimensional_force*(nondimensional_well_parameter*nondimensional_force).cosh() - (nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()) - 3.0 + 3.0*nondimensional_force.ln() - (nondimensional_well_parameter*nondimensional_force*(nondimensional_well_parameter*nondimensional_force).cosh() - (nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()).ln() - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional force, parameterized by the number of links, link length, and well width.
pub fn nondimensional_relative_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, well_width: &f64, nondimensional_force: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, link_length, &1.0, well_width, nondimensional_force, &300.0) - nondimensional_helmholtz_free_energy(number_of_links, link_length, &1.0, well_width, &ZERO, &300.0)
}

/// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional force, parameterized by the link length and well width.
pub fn nondimensional_relative_helmholtz_free_energy_per_link(link_length: &f64, well_width: &f64, nondimensional_force: &f64) -> f64
{
    nondimensional_helmholtz_free_energy_per_link(link_length, &1.0, well_width, nondimensional_force, &300.0) - nondimensional_helmholtz_free_energy_per_link(link_length, &1.0, well_width, &ZERO, &300.0)
}

/// The implemented functionality of the thermodynamics of the SWFJC model in the isotensional ensemble approximated using a Legendre transformation.
impl SWFJC
{
    /// Initializes and returns an instance of the thermodynamics of the SWFJC model in the isotensional ensemble approximated using a Legendre transformation.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, well_width: f64) -> Self
    {
        SWFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            well_width
        }
    }
    /// The Helmholtz free energy as a function of the applied force and temperature.
    pub fn helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.well_width, force, temperature)
    }
    /// The Helmholtz free energy per link as a function of the applied force and temperature.
    pub fn helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy_per_link(&self.link_length, &self.hinge_mass, &self.well_width, force, temperature)
    }
    /// The relative Helmholtz free energy as a function of the applied force and temperature.
    pub fn relative_helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.well_width, force, temperature)
    }
    /// The relative Helmholtz free energy per link as a function of the applied force and temperature.
    pub fn relative_helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy_per_link(&self.link_length, &self.well_width, force, temperature)
    }
    /// The nondimensional Helmholtz free energy as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.well_width, nondimensional_force, temperature)
    }
    /// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy_per_link(&self.link_length, &self.hinge_mass, &self.well_width, nondimensional_force, temperature)
    }
    /// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional force.
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.well_width, nondimensional_force)
    }
    /// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional force.
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy_per_link(&self.link_length, &self.well_width, nondimensional_force)
    }
}
