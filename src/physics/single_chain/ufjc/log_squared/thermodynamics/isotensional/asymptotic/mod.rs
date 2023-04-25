#[cfg(feature = "extern")]
pub mod ex;

#[cfg(feature = "python")]
pub mod py;

mod test;

/// The log-squared link potential freely-jointed chain (log-squared-FJC) model thermodynamics in the isotensional ensemble approximated using a reduced asymptotic approach.
pub mod reduced;

/// The log-squared link potential freely-jointed chain (log-squared-FJC) model thermodynamics in the isotensional ensemble approximated using an asymptotic approach and a Legendre transformation.
pub mod legendre;

use std::f64::consts::PI;
use crate::math::lambert_w;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT,
    single_chain::ZERO
};

/// The structure of the log-squared-FJC model thermodynamics in the isotensional ensemble approximated using an asymptotic approach.
pub struct LOGSQUAREDFJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The stiffness of each link in the chain in units of J/(mol⋅nm^2).
    pub link_stiffness: f64,

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using a reduced asymptotic approach.
    pub reduced: self::reduced::LOGSQUAREDFJC,

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using an asymptotic approach and a Legendre transformation.
    pub legendre: self::legendre::LOGSQUAREDFJC
}

/// The expected end-to-end length as a function of the applied force and temperature, parameterized by the number of links, link length, and link stiffness.
pub fn end_to_end_length(number_of_links: &u8, link_length: &f64, link_stiffness: &f64, force: &f64, temperature: &f64) -> f64
{
    link_length*nondimensional_end_to_end_length(number_of_links, &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(force*link_length/BOLTZMANN_CONSTANT/temperature))
}

/// The expected end-to-end length per link as a function of the applied force and temperature, parameterized by the link length and link stiffness.
pub fn end_to_end_length_per_link(link_length: &f64, link_stiffness: &f64, force: &f64, temperature: &f64) -> f64
{
    link_length*nondimensional_end_to_end_length_per_link(&(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(force*link_length/BOLTZMANN_CONSTANT/temperature))
}

/// The expected nondimensional end-to-end length as a function of the applied nondimensional force, parameterized by the number of links and nondimensional link stiffness.
pub fn nondimensional_end_to_end_length(number_of_links: &u8, nondimensional_link_stiffness: &f64, nondimensional_force: &f64) -> f64
{
    (*number_of_links as f64)*nondimensional_end_to_end_length_per_link(nondimensional_link_stiffness, nondimensional_force)
}

/// The expected nondimensional end-to-end length per link as a function of the applied nondimensional force, parameterized by the nondimensional link stiffness.
pub fn nondimensional_end_to_end_length_per_link(nondimensional_link_stiffness: &f64, nondimensional_force: &f64) -> f64
{
    1.0/nondimensional_force.tanh() - 1.0/nondimensional_force + nondimensional_force/nondimensional_link_stiffness*((nondimensional_force.tanh() - 1.0/nondimensional_force.tanh() + 1.0/nondimensional_force)/(0.4*nondimensional_force.tanh() + nondimensional_force/nondimensional_link_stiffness)) + 1.0/lambert_w(&(-nondimensional_force/nondimensional_link_stiffness)).exp() - 1.0
}

/// The Gibbs free energy as a function of the applied force and temperature, parameterized by the number of links, link length, hinge mass, and link stiffness.
pub fn gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, link_stiffness: &f64, force: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*nondimensional_gibbs_free_energy(number_of_links, link_length, hinge_mass, &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(force*link_length/BOLTZMANN_CONSTANT/temperature), temperature)
}

/// The Gibbs free energy per link as a function of the applied force and temperature, parameterized by the link length, hinge mass, and link stiffness.
pub fn gibbs_free_energy_per_link(link_length: &f64, hinge_mass: &f64, link_stiffness: &f64, force: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*nondimensional_gibbs_free_energy_per_link(link_length, hinge_mass, &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(force*link_length/BOLTZMANN_CONSTANT/temperature), temperature)
}

/// The relative Gibbs free energy as a function of the applied force and temperature, parameterized by the number of links, link length, and link stiffness.
pub fn relative_gibbs_free_energy(number_of_links: &u8, link_length: &f64, link_stiffness: &f64, force: &f64, temperature: &f64) -> f64
{
    gibbs_free_energy(number_of_links, link_length, &1.0, link_stiffness, force, temperature) - gibbs_free_energy(number_of_links, link_length, &1.0, link_stiffness, &(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), temperature)
}

/// The relative Gibbs free energy per link as a function of the applied force and temperature, parameterized by the link length and link stiffness,.
pub fn relative_gibbs_free_energy_per_link(link_length: &f64, link_stiffness: &f64, force: &f64, temperature: &f64) -> f64
{
    gibbs_free_energy_per_link(link_length, &1.0, link_stiffness, force, temperature) - gibbs_free_energy_per_link(link_length, &1.0, link_stiffness, &(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), temperature)
}

/// The nondimensional Gibbs free energy as a function of the applied nondimensional force and temperature, parameterized by the number of links, link length, hinge mass, and nondimensional link stiffness.
pub fn nondimensional_gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_link_stiffness: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    (*number_of_links as f64)*nondimensional_gibbs_free_energy_per_link(link_length, hinge_mass, nondimensional_link_stiffness, nondimensional_force, temperature)
}

/// The nondimensional Gibbs free energy per link as a function of the applied nondimensional force and temperature, parameterized by the link length, hinge mass, and nondimensional link stiffness.
pub fn nondimensional_gibbs_free_energy_per_link(link_length: &f64, hinge_mass: &f64, nondimensional_link_stiffness: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    let lambda = 1.0/lambert_w(&(-nondimensional_force/nondimensional_link_stiffness)).exp();
    -(nondimensional_force.sinh()/nondimensional_force).ln() - (1.0 + nondimensional_force/nondimensional_force.tanh()/nondimensional_link_stiffness/0.4).ln() + 0.5*nondimensional_link_stiffness*lambda.ln().powi(2) - nondimensional_force*(lambda - 1.0) - 0.5*(2.0*PI*link_length.powi(2)/nondimensional_link_stiffness).ln() - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional relative Gibbs free energy as a function of the applied nondimensional force, parameterized by the number of links and nondimensional link stiffness.
pub fn nondimensional_relative_gibbs_free_energy(number_of_links: &u8, nondimensional_link_stiffness: &f64, nondimensional_force: &f64) -> f64
{
    nondimensional_gibbs_free_energy(number_of_links, &1.0, &1.0, nondimensional_link_stiffness, nondimensional_force, &300.0) - nondimensional_gibbs_free_energy(number_of_links, &1.0, &1.0, nondimensional_link_stiffness, &ZERO, &300.0)
}

/// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional force, parameterized by the nondimensional link stiffness.
pub fn nondimensional_relative_gibbs_free_energy_per_link(nondimensional_link_stiffness: &f64, nondimensional_force: &f64) -> f64
{
    nondimensional_gibbs_free_energy_per_link(&1.0, &1.0, nondimensional_link_stiffness, nondimensional_force, &300.0) - nondimensional_gibbs_free_energy_per_link(&1.0, &1.0, nondimensional_link_stiffness, &ZERO, &300.0)
}

/// The implemented functionality of the log-squared-FJC model thermodynamics in the isotensional ensemble approximated using an asymptotic approach.
impl LOGSQUAREDFJC
{
    /// Initializes and returns an instance of the log-squared-FJC model thermodynamics in the isotensional ensemble approximated using an asymptotic approach.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        LOGSQUAREDFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            reduced: self::reduced::LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness),
            legendre: self::legendre::LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness)
        }
    }
    /// The expected end-to-end length as a function of the applied force and temperature.
    pub fn end_to_end_length(&self, force: &f64, temperature: &f64) -> f64
    {
        end_to_end_length(&self.number_of_links, &self.link_length, &self.link_stiffness, force, temperature)
    }
    /// The expected end-to-end length per link as a function of the applied force and temperature.
    pub fn end_to_end_length_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        end_to_end_length_per_link(&self.link_length, &self.link_stiffness, force, temperature)
    }
    /// The expected nondimensional end-to-end length as a function of the applied nondimensional force.
    pub fn nondimensional_end_to_end_length(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_end_to_end_length(&self.number_of_links, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), nondimensional_force)
    }
    /// The expected nondimensional end-to-end length per link as a function of the applied nondimensional force.
    pub fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_end_to_end_length_per_link(&(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), nondimensional_force)
    }
    /// The Gibbs free energy as a function of the applied force and temperature.
    pub fn gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.link_stiffness, force, temperature)
    }
    /// The Gibbs free energy per link as a function of the applied force and temperature.
    pub fn gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        gibbs_free_energy_per_link(&self.link_length, &self.hinge_mass, &self.link_stiffness, force, temperature)
    }
    /// The relative Gibbs free energy as a function of the applied force and temperature.
    pub fn relative_gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.link_stiffness, force, temperature)
    }
    /// The relative Gibbs free energy per link as a function of the applied force and temperature.
    pub fn relative_gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_gibbs_free_energy_per_link(&self.link_length, &self.link_stiffness, force, temperature)
    }
    /// The nondimensional Gibbs free energy as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), nondimensional_force, temperature)
    }
    /// The nondimensional Gibbs free energy per link as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_gibbs_free_energy_per_link(&self.link_length, &self.hinge_mass, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), nondimensional_force, temperature)
    }
    /// The nondimensional relative Gibbs free energy as a function of the applied nondimensional force.
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_relative_gibbs_free_energy(&self.number_of_links, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), nondimensional_force)
    }
    /// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional force.
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_relative_gibbs_free_energy_per_link(&(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), nondimensional_force)
    }
}
