#[cfg(feature = "extern")]
pub mod ex;

#[cfg(feature = "python")]
pub mod py;

mod test;

/// The freely-jointed chain (FJC) model thermodynamics in the modified canonical ensemble approximated using an asymptotic approach.
pub mod asymptotic;

use super::treloar_sum_0_with_prefactor;
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

/// The structure of the thermodynamics of the FJC model in the modified canonical ensemble.
pub struct FJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using an asymptotic approach.
    pub asymptotic: asymptotic::FJC
}

/// The expected end-to-end length as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links and link length.
pub fn end_to_end_length(number_of_links: &u8, link_length: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    potential_distance - force(number_of_links, link_length, potential_distance, potential_stiffness, temperature)/potential_stiffness
}

/// The expected end-to-end length per link as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links and link length.
pub fn end_to_end_length_per_link(number_of_links: &u8, link_length: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    (potential_distance - force(number_of_links, link_length, potential_distance, potential_stiffness, temperature)/potential_stiffness)/(*number_of_links as f64)
}

/// The expected nondimensional end-to-end length as a function of the applied nondimensional potential distance and nondimensional potential stiffness, parameterized by the number of links.
pub fn nondimensional_end_to_end_length(number_of_links: &u8, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
{
    (*number_of_links as f64)*nondimensional_potential_distance - nondimensional_force(number_of_links, nondimensional_potential_distance, nondimensional_potential_stiffness)/nondimensional_potential_stiffness
}

/// The expected nondimensional end-to-end length per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness, parameterized by the number of links.
pub fn nondimensional_end_to_end_length_per_link(number_of_links: &u8, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
{
    nondimensional_potential_distance - nondimensional_force(number_of_links, nondimensional_potential_distance, nondimensional_potential_stiffness)/nondimensional_potential_stiffness/(*number_of_links as f64)
}

/// The expected force as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links and link length.
pub fn force(number_of_links: &u8, link_length: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force(number_of_links, &(*potential_distance/((*number_of_links as f64)*link_length)), &(potential_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature))
}

/// The expected nondimensional force as a function of the applied nondimensional potential distance and nondimensional potential stiffness, parameterized by the number of links.
pub fn nondimensional_force(number_of_links: &u8, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
{
    let number_of_links_f64 = *number_of_links as f64;
    let number_of_links_squared_times_nondimensional_potential_stiffness = number_of_links_f64.powi(2)*nondimensional_potential_stiffness;
    let integrand_numerator = |nondimensional_end_to_end_length_per_link: &f64|
    {
        4.0*PI*nondimensional_end_to_end_length_per_link.powi(2)*treloar_sum_0_with_prefactor(number_of_links, nondimensional_end_to_end_length_per_link)*((number_of_links_squared_times_nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link)*(-0.5*number_of_links_squared_times_nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link).powi(2)).exp() - number_of_links_squared_times_nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link)*(-0.5*number_of_links_squared_times_nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link).powi(2)).exp())/(2.0*number_of_links_squared_times_nondimensional_potential_stiffness*nondimensional_potential_distance*nondimensional_end_to_end_length_per_link) + ((-0.5*number_of_links_squared_times_nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link).powi(2)).exp() - (-0.5*number_of_links_squared_times_nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link).powi(2)).exp())/(2.0*number_of_links_squared_times_nondimensional_potential_stiffness*nondimensional_potential_distance.powi(2)*nondimensional_end_to_end_length_per_link))
    };
    let integrand_denominator = |nondimensional_end_to_end_length_per_link: &f64|
    {
        4.0*PI*nondimensional_end_to_end_length_per_link.powi(2)*treloar_sum_0_with_prefactor(number_of_links, nondimensional_end_to_end_length_per_link)*((-0.5*number_of_links_squared_times_nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link).powi(2)).exp() - (-0.5*number_of_links_squared_times_nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link).powi(2)).exp())/(2.0*number_of_links_squared_times_nondimensional_potential_stiffness*nondimensional_potential_distance*nondimensional_end_to_end_length_per_link)
    };
    integrate_1d(&integrand_numerator, &ZERO, &ONE, &POINTS)/integrate_1d(&integrand_denominator, &ZERO, &ONE, &POINTS)/number_of_links_f64
}

/// The Helmholtz free energy as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, &(potential_distance/((*number_of_links as f64)*link_length)), &(potential_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), temperature)
}

/// The Helmholtz free energy per link as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*nondimensional_helmholtz_free_energy_per_link(number_of_links, link_length, hinge_mass, &(*potential_distance/((*number_of_links as f64)*link_length)), &(potential_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), temperature)
}

/// The relative Helmholtz free energy as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links and link length.
pub fn relative_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*nondimensional_relative_helmholtz_free_energy(number_of_links, &(*potential_distance/((*number_of_links as f64)*link_length)), &(potential_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature))
}

/// The relative Helmholtz free energy per link as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links and link length.
pub fn relative_helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*nondimensional_relative_helmholtz_free_energy_per_link(number_of_links, &(*potential_distance/((*number_of_links as f64)*link_length)), &(potential_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature))
}

/// The nondimensional Helmholtz free energy as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
{
    let number_of_links_f64 = *number_of_links as f64;
    let number_of_links_squared_times_nondimensional_potential_stiffness = number_of_links_f64.powi(2)*nondimensional_potential_stiffness;
    let integrand = |nondimensional_end_to_end_length_per_link: &f64|
    {
        4.0*PI*nondimensional_end_to_end_length_per_link.powi(2)*treloar_sum_0_with_prefactor(number_of_links, nondimensional_end_to_end_length_per_link)*((-0.5*number_of_links_squared_times_nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link).powi(2)).exp() - (-0.5*number_of_links_squared_times_nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link).powi(2)).exp())/(2.0*number_of_links_squared_times_nondimensional_potential_stiffness*nondimensional_potential_distance*nondimensional_end_to_end_length_per_link)
    };
    let nondimensional_configurational_partition_function = integrate_1d(&integrand, &ZERO, &ONE, &POINTS);
    -nondimensional_configurational_partition_function.ln() - (number_of_links_f64 - 1.0)*(8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, nondimensional_potential_distance, nondimensional_potential_stiffness, temperature)/(*number_of_links as f64)
}

/// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional potential distance and nondimensional potential stiffness, parameterized by the number of links.
pub fn nondimensional_relative_helmholtz_free_energy(number_of_links: &u8, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, &1.0, &1.0, nondimensional_potential_distance, nondimensional_potential_stiffness, &300.0) - nondimensional_helmholtz_free_energy(number_of_links, &1.0, &1.0, &ZERO, nondimensional_potential_stiffness, &300.0)
}

/// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness, parameterized by the number of links.
pub fn nondimensional_relative_helmholtz_free_energy_per_link(number_of_links: &u8, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy(number_of_links, nondimensional_potential_distance, nondimensional_potential_stiffness)/(*number_of_links as f64)
}

/// The Gibbs free energy as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    helmholtz_free_energy(number_of_links, link_length, hinge_mass, potential_distance, potential_stiffness, temperature) - 0.5*potential_stiffness*potential_distance.powi(2)
}

/// The Gibbs free energy epr link as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    helmholtz_free_energy_per_link(number_of_links, link_length, hinge_mass, potential_distance, potential_stiffness, temperature) - 0.5*potential_stiffness*potential_distance.powi(2)/(*number_of_links as f64)
}

/// The relative Gibbs free energy as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links and link length.
pub fn relative_gibbs_free_energy(number_of_links: &u8, link_length: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    relative_helmholtz_free_energy(number_of_links, link_length, potential_distance, potential_stiffness, temperature) - 0.5*potential_stiffness*potential_distance.powi(2)
}

/// The relative Gibbs free energy per link as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links and link length.
pub fn relative_gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    relative_helmholtz_free_energy_per_link(number_of_links, link_length, potential_distance, potential_stiffness, temperature) - 0.5*potential_stiffness*potential_distance.powi(2)/(*number_of_links as f64)
}

/// The nondimensional Gibbs free energy as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, nondimensional_potential_distance, nondimensional_potential_stiffness, temperature) - 0.5*(*number_of_links as f64).powi(2)*nondimensional_potential_stiffness*nondimensional_potential_distance.powi(2)
}

/// The nondimensional Gibbs free energy per link as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy_per_link(number_of_links, link_length, hinge_mass, nondimensional_potential_distance, nondimensional_potential_stiffness, temperature) - 0.5*(*number_of_links as f64)*nondimensional_potential_stiffness*nondimensional_potential_distance.powi(2)
}

/// The nondimensional relative Gibbs free energy as a function of the applied nondimensional potential distance and nondimensional potential stiffness, parameterized by the number of links.
pub fn nondimensional_relative_gibbs_free_energy(number_of_links: &u8, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy(number_of_links, nondimensional_potential_distance, nondimensional_potential_stiffness) - 0.5*(*number_of_links as f64).powi(2)*nondimensional_potential_stiffness*nondimensional_potential_distance.powi(2)
}

/// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness, parameterized by the number of links.
pub fn nondimensional_relative_gibbs_free_energy_per_link(number_of_links: &u8, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy_per_link(number_of_links, nondimensional_potential_distance, nondimensional_potential_stiffness) - 0.5*(*number_of_links as f64)*nondimensional_potential_stiffness*nondimensional_potential_distance.powi(2)
}

/// The implemented functionality of the thermodynamics of the FJC model in the modified canonical ensemble.
impl FJC
{
    /// Initializes and returns an instance of the thermodynamics of the FJC model in the modified canonical ensemble.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            asymptotic: asymptotic::FJC::init(number_of_links, link_length, hinge_mass)
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
    /// The expected force as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn force(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        force(&self.number_of_links, &self.link_length, potential_distance, potential_stiffness, temperature)
    }
    /// The expected nondimensional force as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    pub fn nondimensional_force(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        nondimensional_force(&self.number_of_links, nondimensional_potential_distance, nondimensional_potential_stiffness)
    }
    /// The Helmholtz free energy as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn helmholtz_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, potential_distance, potential_stiffness, temperature)
    }
    /// The Helmholtz free energy per link as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn helmholtz_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, potential_distance, potential_stiffness, temperature)
    }
    /// The relative Helmholtz free energy as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn relative_helmholtz_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy(&self.number_of_links, &self.link_length, potential_distance, potential_stiffness, temperature)
    }
    /// The relative Helmholtz free energy per link as a function of the applied potential distance, potential stiffness, and temperature.
    pub fn relative_helmholtz_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, potential_distance, potential_stiffness, temperature)
    }
    /// The nondimensional Helmholtz free energy as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, nondimensional_potential_distance, nondimensional_potential_stiffness, temperature)
    }
    /// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature.
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, nondimensional_potential_distance, nondimensional_potential_stiffness, temperature)
    }
    /// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy(&self.number_of_links, nondimensional_potential_distance, nondimensional_potential_stiffness)
    }
    /// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional potential distance and nondimensional potential stiffness.
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy_per_link(&self.number_of_links, nondimensional_potential_distance, nondimensional_potential_stiffness)
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
