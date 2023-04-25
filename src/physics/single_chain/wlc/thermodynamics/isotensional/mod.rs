#[cfg(feature = "extern")]
pub mod ex;

#[cfg(feature = "python")]
pub mod py;

mod test;

/// The worm-like chain (WLC) model thermodynamics in the isotensional ensemble approximated using a Legendre transformation.
pub mod legendre;

use super::isometric::nondimensional_helmholtz_free_energy;

use std::f64::consts::PI;
use crate::math::integrate_1d;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::
{
    ONE,
    ZERO,
    POINTS
};

/// The structure of the thermodynamics of the WLC model in the isotensional ensemble.
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

    nondimensional_persistance_length: f64,

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using a Legendre transformation.
    pub legendre: self::legendre::WLC
}

/// The expected end-to-end length as a function of the applied force and temperature, parameterized by the number of links, link length, and persistance length.
pub fn end_to_end_length(number_of_links: &u8, link_length: &f64, persistance_length: &f64, force: &f64, temperature: &f64) -> f64
{
    link_length*nondimensional_end_to_end_length(number_of_links, &(persistance_length/((*number_of_links as f64)*link_length)), &(force*link_length/BOLTZMANN_CONSTANT/temperature))
}

/// The expected end-to-end length per link as a function of the applied force and temperature, parameterized by the number of links, link length, and persistance length.
pub fn end_to_end_length_per_link(number_of_links: &u8, link_length: &f64, persistance_length: &f64, force: &f64, temperature: &f64) -> f64
{
    link_length*nondimensional_end_to_end_length_per_link(number_of_links, &(persistance_length/((*number_of_links as f64)*link_length)), &(force*link_length/BOLTZMANN_CONSTANT/temperature))
}

/// The expected nondimensional end-to-end length as a function of the applied nondimensional force, parameterized by the number of links and nondimensional persistance length.
pub fn nondimensional_end_to_end_length(number_of_links: &u8, nondimensional_persistance_length: &f64, nondimensional_force: &f64) -> f64
{
    let integrand_numerator = |nondimensional_end_to_end_length_per_link: &f64| (-nondimensional_relative_gibbs_free_energy(nondimensional_persistance_length, nondimensional_end_to_end_length_per_link)).exp()*((*number_of_links as f64)*nondimensional_end_to_end_length_per_link*((*number_of_links as f64)*nondimensional_force*nondimensional_end_to_end_length_per_link).cosh()*nondimensional_end_to_end_length_per_link - ((*number_of_links as f64)*nondimensional_force*nondimensional_end_to_end_length_per_link).sinh()*nondimensional_end_to_end_length_per_link/nondimensional_force)/nondimensional_force;
    let integrand_denominator = |nondimensional_end_to_end_length_per_link: &f64| (-nondimensional_relative_gibbs_free_energy(nondimensional_persistance_length, nondimensional_end_to_end_length_per_link)).exp()*((*number_of_links as f64)*nondimensional_force*nondimensional_end_to_end_length_per_link).sinh()*nondimensional_end_to_end_length_per_link/nondimensional_force;
    integrate_1d(&integrand_numerator, &ZERO, &ONE, &POINTS).ln()/integrate_1d(&integrand_denominator, &ZERO, &ONE, &POINTS).ln()
}

/// The expected nondimensional end-to-end length per link as a function of the applied nondimensional force, parameterized by the number of links and nondimensional persistance length.
pub fn nondimensional_end_to_end_length_per_link(number_of_links: &u8, nondimensional_persistance_length: &f64, nondimensional_force: &f64) -> f64
{
    (*number_of_links as f64)*nondimensional_end_to_end_length(number_of_links, nondimensional_persistance_length, nondimensional_force)
}

/// The Gibbs free energy as a function of the applied force and temperature, parameterized by the number of links, link length, hinge mass, and persistance length.
pub fn gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, persistance_length: &f64, force: &f64, temperature: &f64) -> f64
{
    nondimensional_gibbs_free_energy(number_of_links, link_length, hinge_mass, &(persistance_length/((*number_of_links as f64)*link_length)), &(force*link_length/BOLTZMANN_CONSTANT/temperature), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The Gibbs free energy per link as a function of the applied force and temperature, parameterized by the number of links, link length, hinge mass, and persistance length.
pub fn gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, persistance_length: &f64, force: &f64, temperature: &f64) -> f64
{
    nondimensional_gibbs_free_energy_per_link(number_of_links, link_length, hinge_mass, &(persistance_length/((*number_of_links as f64)*link_length)), &(force*link_length/BOLTZMANN_CONSTANT/temperature), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The relative Gibbs free energy as a function of the applied force and temperature, parameterized by the number of links, link length, and persistance length.
pub fn relative_gibbs_free_energy(number_of_links: &u8, link_length: &f64, persistance_length: &f64, force: &f64, temperature: &f64) -> f64
{
    nondimensional_relative_gibbs_free_energy(&(persistance_length/((*number_of_links as f64)*link_length)), &(force*link_length/BOLTZMANN_CONSTANT/temperature))*BOLTZMANN_CONSTANT*temperature
}

/// The relative Gibbs free energy per link as a function of the applied force and temperature, parameterized by the number of links, link length, and persistance length.
pub fn relative_gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, persistance_length: &f64, force: &f64, temperature: &f64) -> f64
{
    nondimensional_relative_gibbs_free_energy_per_link(&(persistance_length/((*number_of_links as f64)*link_length)), &(force*link_length/BOLTZMANN_CONSTANT/temperature))*BOLTZMANN_CONSTANT*temperature
}

/// The nondimensional Gibbs free energy as a function of the applied nondimensional force and temperature, parameterized by the number of links, link length, hinge mass, and nondimensional persistance length.
pub fn nondimensional_gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_persistance_length: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    let integrand = |nondimensional_end_to_end_length_per_link: &f64|
    {
        (-nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, nondimensional_persistance_length, nondimensional_end_to_end_length_per_link, temperature)).exp()*((*number_of_links as f64)*nondimensional_force*nondimensional_end_to_end_length_per_link).sinh()*(*number_of_links as f64)*nondimensional_end_to_end_length_per_link/nondimensional_force
    };
    (4.0*PI*integrate_1d(&integrand, &ZERO, &ONE, &POINTS)).ln() - (4.0*(-1.0/nondimensional_persistance_length).exp().acos().sin()*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional Gibbs free energy per link as a function of the applied nondimensional force and temperature, parameterized by the number of links, link length, and hinge mass, and nondimensional persistance length.
pub fn nondimensional_gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_persistance_length: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    nondimensional_gibbs_free_energy(number_of_links, link_length, hinge_mass, nondimensional_persistance_length, nondimensional_force, temperature)/(*number_of_links as f64)
}

/// The nondimensional relative Gibbs free energy as a function of the applied nondimensional force, parameterized by the nondimensional persistance length.
pub fn nondimensional_relative_gibbs_free_energy(nondimensional_persistance_length: &f64, nondimensional_force: &f64) -> f64
{
    nondimensional_gibbs_free_energy(&8, &1.0, &1.0, nondimensional_persistance_length, nondimensional_force, &300.0) - nondimensional_gibbs_free_energy(&8, &1.0, &1.0, nondimensional_persistance_length, &ZERO, &300.0)
}

/// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional force, parameterized by the nondimensional persistance length.
pub fn nondimensional_relative_gibbs_free_energy_per_link(nondimensional_persistance_length: &f64, nondimensional_force: &f64) -> f64
{
    nondimensional_gibbs_free_energy_per_link(&8, &1.0, &1.0, nondimensional_persistance_length, nondimensional_force, &300.0) - nondimensional_gibbs_free_energy_per_link(&8, &1.0, &1.0, nondimensional_persistance_length, &ZERO, &300.0)
}

/// The implemented functionality of the thermodynamics of the WLC model in the isotensional ensemble.
impl WLC
{
    /// Initializes and returns an instance of the thermodynamics of the WLC model in the isotensional ensemble.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, persistance_length: f64) -> Self
    {
        WLC
        {
            hinge_mass,
            link_length,
            number_of_links,
            persistance_length,
            nondimensional_persistance_length: persistance_length/(number_of_links as f64)/link_length,
            legendre: self::legendre::WLC::init(number_of_links, link_length, hinge_mass, persistance_length),
        }
    }
    /// The expected end-to-end length as a function of the applied force and temperature.
    pub fn end_to_end_length(&self, force: &f64, temperature: &f64) -> f64
    {
        end_to_end_length(&self.number_of_links, &self.link_length, &self.persistance_length, force, temperature)
    }
    /// The expected end-to-end length per link as a function of the applied force and temperature.
    pub fn end_to_end_length_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        end_to_end_length_per_link(&self.number_of_links, &self.link_length, &self.persistance_length, force, temperature)
    }
    /// The expected nondimensional end-to-end length as a function of the applied nondimensional force.
    pub fn nondimensional_end_to_end_length(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_end_to_end_length(&self.number_of_links, &self.nondimensional_persistance_length, nondimensional_force)
    }
    /// The expected nondimensional end-to-end length per link as a function of the applied nondimensional force.
    pub fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_end_to_end_length_per_link(&self.number_of_links, &self.nondimensional_persistance_length, nondimensional_force)
    }
    /// The Gibbs free energy as a function of the applied force and temperature.
    pub fn gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.persistance_length, force, temperature)
    }
    /// The Gibbs free energy per link as a function of the applied force and temperature.
    pub fn gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.persistance_length, force, temperature)
    }
    /// The relative Gibbs free energy as a function of the applied force and temperature.
    pub fn relative_gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.persistance_length, force, temperature)
    }
    /// The relative Gibbs free energy per link as a function of the applied force and temperature.
    pub fn relative_gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.persistance_length, force, temperature)
    }
    /// The nondimensional Gibbs free energy as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.nondimensional_persistance_length, nondimensional_force, temperature)
    }
    /// The nondimensional Gibbs free energy per link as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.nondimensional_persistance_length, nondimensional_force, temperature)
    }
    /// The nondimensional relative Gibbs free energy as a function of the applied nondimensional force.
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_relative_gibbs_free_energy(&self.nondimensional_persistance_length, nondimensional_force)
    }
    /// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional force.
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: &f64) -> f64
    {
        nondimensional_relative_gibbs_free_energy_per_link(&self.nondimensional_persistance_length, nondimensional_force)
    }
}