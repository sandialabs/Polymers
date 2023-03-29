#[cfg(feature = "python")]
pub mod py;

mod test;

/// The Morse link potential freely-jointed chain (Morse-FJC) model thermodynamics in the isotensional ensemble approximated using a reduced asymptotic approach and a Legendre transformation.
pub mod legendre;

use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT,
    single_chain::ZERO
};

/// The structure of the Morse-FJC model thermodynamics in the isotensional ensemble approximated using a reduced asymptotic approach.
pub struct MORSEFJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The stiffness of each link in the chain in units of J/(mol⋅nm^2).
    pub link_stiffness: f64,

    /// The energy of each link in the chain in units of J/mol.
    pub link_energy: f64,

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using a reduced asymptotic approach and a Legendre transformation.
    pub legendre: self::legendre::MORSEFJC
}

/// The expected end-to-end length as a function of the applied force and temperature, parameterized by the number of links, link length, link stiffness, and link energy.
pub fn end_to_end_length(number_of_links: &u8, link_length: &f64, link_stiffness: &f64, link_energy: &f64, force: &f64, temperature: &f64) -> f64
{
    link_length*nondimensional_end_to_end_length(number_of_links, &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(force*link_length/BOLTZMANN_CONSTANT/temperature))
}

/// The expected end-to-end length per link as a function of the applied force and temperature, parameterized by the link length, link stiffness, and link energy.
pub fn end_to_end_length_per_link(link_length: &f64, link_stiffness: &f64, link_energy: &f64, force: &f64, temperature: &f64) -> f64
{
    link_length*nondimensional_end_to_end_length_per_link(&(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(force*link_length/BOLTZMANN_CONSTANT/temperature))
}

/// The expected nondimensional end-to-end length as a function of the applied nondimensional force, parameterized by the number of links, the nondimensional link stiffness, and nondimensional link energy.
pub fn nondimensional_end_to_end_length(number_of_links: &u8, nondimensional_link_stiffness: &f64, nondimensional_link_energy: &f64, nondimensional_force: &f64) -> f64
{
    (*number_of_links as f64)*nondimensional_end_to_end_length_per_link(nondimensional_link_stiffness, nondimensional_link_energy, nondimensional_force)
}

/// The expected nondimensional end-to-end length per link as a function of the applied nondimensional force, parameterized by the nondimensional link stiffness and nondimensional link energy.
pub fn nondimensional_end_to_end_length_per_link(nondimensional_link_stiffness: &f64, nondimensional_link_energy: &f64, nondimensional_force: &f64) -> f64
{
    1.0/nondimensional_force.tanh() - 1.0/nondimensional_force + (2.0*nondimensional_link_energy/nondimensional_link_stiffness).sqrt()*(2.0/(1.0 + (1.0 - nondimensional_force/(nondimensional_link_energy*nondimensional_link_stiffness/8.0).sqrt()).sqrt())).ln()
}

/// The Gibbs free energy as a function of the applied force and temperature, parameterized by the number of links, link length, hinge mass, link stiffness, and link energy.
pub fn gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, link_stiffness: &f64, link_energy: &f64, force: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*nondimensional_gibbs_free_energy(number_of_links, link_length, hinge_mass, &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(force*link_length/BOLTZMANN_CONSTANT/temperature), temperature)
}

/// The Gibbs free energy per link as a function of the applied force and temperature, parameterized by the link length, hinge mass, link stiffness, and link energy.
pub fn gibbs_free_energy_per_link(link_length: &f64, hinge_mass: &f64, link_stiffness: &f64, link_energy: &f64, force: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*nondimensional_gibbs_free_energy_per_link(link_length, hinge_mass, &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(force*link_length/BOLTZMANN_CONSTANT/temperature), temperature)
}

/// The relative Gibbs free energy as a function of the applied force and temperature, parameterized by the number of links, link length, link stiffness, and link energy.
pub fn relative_gibbs_free_energy(number_of_links: &u8, link_length: &f64, link_stiffness: &f64, link_energy: &f64, force: &f64, temperature: &f64) -> f64
{
    gibbs_free_energy(number_of_links, link_length, &1.0, link_stiffness, link_energy, force, temperature) - gibbs_free_energy(number_of_links, link_length, &1.0, link_stiffness, link_energy, &(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), temperature)
}

/// The relative Gibbs free energy per link as a function of the applied force and temperature, parameterized by the link length, link stiffness, and link energy.
pub fn relative_gibbs_free_energy_per_link(link_length: &f64, link_stiffness: &f64, link_energy: &f64, force: &f64, temperature: &f64) -> f64
{
    gibbs_free_energy_per_link(link_length, &1.0, link_stiffness, link_energy, force, temperature) - gibbs_free_energy_per_link(link_length, &1.0, link_stiffness, link_energy, &(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), temperature)
}

/// The nondimensional Gibbs free energy as a function of the applied nondimensional force and temperature, parameterized by the number of links, link length, hinge mass, nondimensional link stiffness, and nondimensional link energy.
pub fn nondimensional_gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_link_stiffness: &f64, nondimensional_link_energy: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    (*number_of_links as f64)*nondimensional_gibbs_free_energy_per_link(link_length, hinge_mass, nondimensional_link_stiffness, nondimensional_link_energy, nondimensional_force, temperature)
}

/// The nondimensional Gibbs free energy per link as a function of the applied nondimensional force and temperature, parameterized by the link length, hinge mass, nondimensional link stiffness, and nondimensional link energy.
pub fn nondimensional_gibbs_free_energy_per_link(link_length: &f64, hinge_mass: &f64, nondimensional_link_stiffness: &f64, nondimensional_link_energy: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    let nondimensional_force_max = (nondimensional_link_energy*nondimensional_link_stiffness/8.0).sqrt();
    -(nondimensional_force.sinh()/nondimensional_force).ln() + nondimensional_link_energy*(1.0 - 0.5*(1.0 + (1.0 - nondimensional_force/nondimensional_force_max).sqrt())).powi(2) - nondimensional_force*(2.0*nondimensional_link_energy/nondimensional_link_stiffness).sqrt()*(2.0/(1.0 + (1.0 - nondimensional_force/nondimensional_force_max).sqrt())).ln() - 0.5*(2.0*PI*link_length.powi(2)/nondimensional_link_stiffness).ln() - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional relative Gibbs free energy as a function of the applied nondimensional force, parameterized by the number of links, nondimensional link stiffness, and nondimensional link energy.
pub fn nondimensional_relative_gibbs_free_energy(number_of_links: &u8, nondimensional_link_stiffness: &f64, nondimensional_link_energy: &f64, nondimensional_force: &f64) -> f64
{
    nondimensional_gibbs_free_energy(number_of_links, &1.0, &1.0, nondimensional_link_stiffness, nondimensional_link_energy, nondimensional_force, &300.0) - nondimensional_gibbs_free_energy(number_of_links, &1.0, &1.0, nondimensional_link_stiffness, nondimensional_link_energy, &ZERO, &300.0)
}

/// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional force, parameterized by the nondimensional link stiffness and nondimensional link energy.
pub fn nondimensional_relative_gibbs_free_energy_per_link(nondimensional_link_stiffness: &f64, nondimensional_link_energy: &f64, nondimensional_force: &f64) -> f64
{
    nondimensional_gibbs_free_energy_per_link(&1.0, &1.0, nondimensional_link_stiffness, nondimensional_link_energy, nondimensional_force, &300.0) - nondimensional_gibbs_free_energy_per_link(&1.0, &1.0, nondimensional_link_stiffness, nondimensional_link_energy, &ZERO, &300.0)
}

/// The implemented functionality of the Morse-FJC model thermodynamics in the isotensional ensemble approximated using a reduced asymptotic approach.
impl MORSEFJC
{
    /// Initializes and returns an instance of the Morse-FJC model thermodynamics in the isotensional ensemble approximated using a reduced asymptotic approach.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64, link_energy: f64) -> Self
    {
        MORSEFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            link_energy,
            legendre: self::legendre::MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy)
        }
    }
    /// The expected end-to-end length as a function of the applied force and temperature.
    pub fn end_to_end_length(&self, force: &f64, temperature: &f64) -> f64
    {
        end_to_end_length(&self.number_of_links, &self.link_length, &self.link_stiffness, &self.link_energy, force, temperature)
    }
    /// The expected end-to-end length per link as a function of the applied force and temperature.
    pub fn end_to_end_length_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        end_to_end_length_per_link(&self.link_length, &self.link_stiffness, &self.link_energy, force, temperature)
    }
    /// The expected nondimensional end-to-end length as a function of the applied nondimensional force.
    pub fn nondimensional_end_to_end_length(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_end_to_end_length(&self.number_of_links, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(self.link_energy/BOLTZMANN_CONSTANT/temperature), nondimensional_force)
    }
    /// The expected nondimensional end-to-end length per link as a function of the applied nondimensional force.
    pub fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_end_to_end_length_per_link(&(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(self.link_energy/BOLTZMANN_CONSTANT/temperature), nondimensional_force)
    }
    /// The Gibbs free energy as a function of the applied force and temperature.
    pub fn gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.link_stiffness, &self.link_energy, force, temperature)
    }
    /// The Gibbs free energy per link as a function of the applied force and temperature.
    pub fn gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        gibbs_free_energy_per_link(&self.link_length, &self.hinge_mass, &self.link_stiffness, &self.link_energy, force, temperature)
    }
    /// The relative Gibbs free energy as a function of the applied force and temperature.
    pub fn relative_gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.link_stiffness, &self.link_energy, force, temperature)
    }
    /// The relative Gibbs free energy per link as a function of the applied force and temperature.
    pub fn relative_gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_gibbs_free_energy_per_link(&self.link_length, &self.link_stiffness, &self.link_energy, force, temperature)
    }
    /// The nondimensional Gibbs free energy as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(self.link_energy/BOLTZMANN_CONSTANT/temperature), nondimensional_force, temperature)
    }
    /// The nondimensional Gibbs free energy per link as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_gibbs_free_energy_per_link(&self.link_length, &self.hinge_mass, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(self.link_energy/BOLTZMANN_CONSTANT/temperature), nondimensional_force, temperature)
    }
    /// The nondimensional relative Gibbs free energy as a function of the applied nondimensional force.
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_relative_gibbs_free_energy(&self.number_of_links, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(self.link_energy/BOLTZMANN_CONSTANT/temperature), nondimensional_force)
    }
    /// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional force.
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_relative_gibbs_free_energy_per_link(&(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(self.link_energy/BOLTZMANN_CONSTANT/temperature), nondimensional_force)
    }
}
