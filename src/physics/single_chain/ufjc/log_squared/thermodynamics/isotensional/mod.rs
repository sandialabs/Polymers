#[cfg(feature = "extern")]
pub mod ex;

#[cfg(feature = "python")]
pub mod py;

mod test;

/// The log-squared link potential freely-jointed chain (log-squared-FJC) model thermodynamics in the isotensional ensemble approximated using an asymptotic approach.
pub mod asymptotic;

/// The log-squared link potential freely-jointed chain (log-squared-FJC) model thermodynamics in the isotensional ensemble approximated using a Legendre transformation.
pub mod legendre;

use std::f64::consts::PI;
use crate::math::integrate_1d;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT,
    single_chain::
    {
        ZERO,
        POINTS
    }
};

/// The structure of the log-squared-FJC model thermodynamics in the isotensional ensemble.
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

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using an asymptotic approach.
    pub asymptotic: self::asymptotic::LOGSQUAREDFJC,

    /// The thermodynamic functions of the model in the isotensional ensemble approximated using a Legendre transformation.
    pub legendre: self::legendre::LOGSQUAREDFJC
}

/// The expected end-to-end length as a function of the applied force and temperature, parameterized by the number of links, link length, link stiffness, and link energy.
pub fn end_to_end_length(number_of_links: &u8, link_length: &f64, link_stiffness: &f64, force: &f64, temperature: &f64) -> f64
{
    link_length*nondimensional_end_to_end_length(number_of_links, &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(force*link_length/BOLTZMANN_CONSTANT/temperature))
}

/// The expected end-to-end length per link as a function of the applied force and temperature, parameterized by the link length, link stiffness, and link energy.
pub fn end_to_end_length_per_link(link_length: &f64, link_stiffness: &f64, force: &f64, temperature: &f64) -> f64
{
    link_length*nondimensional_end_to_end_length_per_link(&(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(force*link_length/BOLTZMANN_CONSTANT/temperature))
}

/// The expected nondimensional end-to-end length as a function of the applied nondimensional force, parameterized by the number of links, the nondimensional link stiffness, and nondimensional link energy.
pub fn nondimensional_end_to_end_length(number_of_links: &u8, nondimensional_link_stiffness: &f64, nondimensional_force: &f64) -> f64
{
    (*number_of_links as f64)*nondimensional_end_to_end_length_per_link(nondimensional_link_stiffness, nondimensional_force)
}

/// The expected nondimensional end-to-end length per link as a function of the applied nondimensional force, parameterized by the nondimensional link stiffness and nondimensional link energy.
pub fn nondimensional_end_to_end_length_per_link(nondimensional_link_stiffness: &f64, nondimensional_force: &f64) -> f64
{
    let nondimensional_link_stretch_max = 1.0_f64.exp();
    let rescaled_partition_function_integrand = |nondimensional_link_stretch: &f64|
    {
        let exponent_1 = nondimensional_force*nondimensional_link_stretch - 0.5*nondimensional_link_stiffness*nondimensional_link_stretch.ln().powi(2) + nondimensional_link_stretch.ln() - nondimensional_force.ln();
        let exponent_2 = exponent_1 - 2.0*nondimensional_force*nondimensional_link_stretch;
        exponent_1.exp() - exponent_2.exp()
    };
    let rescaled_partition_function = integrate_1d(&rescaled_partition_function_integrand, &ZERO, &nondimensional_link_stretch_max, &POINTS);
    let nondimensional_end_to_end_length_per_link_integrand = |nondimensional_link_stretch: &f64|
    {
        let exponent_1 = nondimensional_force*nondimensional_link_stretch - 0.5*nondimensional_link_stiffness*nondimensional_link_stretch.ln().powi(2) + 2.0*nondimensional_link_stretch.ln() - nondimensional_force.ln();
        let exponent_2 = exponent_1 - 2.0*nondimensional_force*nondimensional_link_stretch;
        let exponent_3 = exponent_1 - nondimensional_link_stretch.ln() - nondimensional_force.ln();
        let exponent_4 = exponent_3 - 2.0*nondimensional_force*nondimensional_link_stretch;
        (exponent_1.exp() + exponent_2.exp() - exponent_3.exp() + exponent_4.exp())/rescaled_partition_function
    };
    integrate_1d(&nondimensional_end_to_end_length_per_link_integrand, &ZERO, &nondimensional_link_stretch_max, &POINTS)
}

/// The Gibbs free energy as a function of the applied force and temperature, parameterized by the number of links, link length, hinge mass, link stiffness, and link energy.
pub fn gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, link_stiffness: &f64, force: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*nondimensional_gibbs_free_energy(number_of_links, link_length, hinge_mass, &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(force*link_length/BOLTZMANN_CONSTANT/temperature), temperature)
}

/// The Gibbs free energy per link as a function of the applied force and temperature, parameterized by the link length, hinge mass, link stiffness, and link energy.
pub fn gibbs_free_energy_per_link(link_length: &f64, hinge_mass: &f64, link_stiffness: &f64, force: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*nondimensional_gibbs_free_energy_per_link(link_length, hinge_mass, &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(force*link_length/BOLTZMANN_CONSTANT/temperature), temperature)
}

/// The relative Gibbs free energy as a function of the applied force and temperature, parameterized by the number of links, link length, link stiffness, and link energy.
pub fn relative_gibbs_free_energy(number_of_links: &u8, link_length: &f64, link_stiffness: &f64, force: &f64, temperature: &f64) -> f64
{
    gibbs_free_energy(number_of_links, link_length, &1.0, link_stiffness, force, temperature) - gibbs_free_energy(number_of_links, link_length, &1.0, link_stiffness, &(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), temperature)
}

/// The relative Gibbs free energy per link as a function of the applied force and temperature, parameterized by the link length, link stiffness, and link energy.
pub fn relative_gibbs_free_energy_per_link(link_length: &f64, link_stiffness: &f64, force: &f64, temperature: &f64) -> f64
{
    gibbs_free_energy_per_link(link_length, &1.0, link_stiffness, force, temperature) - gibbs_free_energy_per_link(link_length, &1.0, link_stiffness, &(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), temperature)
}

/// The nondimensional Gibbs free energy as a function of the applied nondimensional force and temperature, parameterized by the number of links, link length, hinge mass, nondimensional link stiffness, and nondimensional link energy.
pub fn nondimensional_gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_link_stiffness: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    (*number_of_links as f64)*nondimensional_gibbs_free_energy_per_link(link_length, hinge_mass, nondimensional_link_stiffness, nondimensional_force, temperature)
}

/// The nondimensional Gibbs free energy per link as a function of the applied nondimensional force and temperature, parameterized by the link length, hinge mass, nondimensional link stiffness, and nondimensional link energy.
pub fn nondimensional_gibbs_free_energy_per_link(link_length: &f64, hinge_mass: &f64, nondimensional_link_stiffness: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    let nondimensional_link_stretch_max = 1.0_f64.exp();
    let rescaled_partition_function_integrand = |nondimensional_link_stretch: &f64|
    {
        let exponent_1 = nondimensional_force*nondimensional_link_stretch - 0.5*nondimensional_link_stiffness*nondimensional_link_stretch.ln().powi(2) + nondimensional_link_stretch.ln() - nondimensional_force.ln();
        let exponent_2 = exponent_1 - 2.0*nondimensional_force*nondimensional_link_stretch;
        exponent_1.exp() - exponent_2.exp()
    };
    let rescaled_partition_function = integrate_1d(&rescaled_partition_function_integrand, &ZERO, &nondimensional_link_stretch_max, &POINTS);
    -rescaled_partition_function.ln() - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional relative Gibbs free energy as a function of the applied nondimensional force, parameterized by the number of links, nondimensional link stiffness, and nondimensional link energy.
pub fn nondimensional_relative_gibbs_free_energy(number_of_links: &u8, nondimensional_link_stiffness: &f64, nondimensional_force: &f64) -> f64
{
    nondimensional_gibbs_free_energy(number_of_links, &1.0, &1.0, nondimensional_link_stiffness, nondimensional_force, &300.0) - nondimensional_gibbs_free_energy(number_of_links, &1.0, &1.0, nondimensional_link_stiffness, &ZERO, &300.0)
}

/// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional force, parameterized by the nondimensional link stiffness and nondimensional link energy.
pub fn nondimensional_relative_gibbs_free_energy_per_link(nondimensional_link_stiffness: &f64, nondimensional_force: &f64) -> f64
{
    nondimensional_gibbs_free_energy_per_link(&1.0, &1.0, nondimensional_link_stiffness, nondimensional_force, &300.0) - nondimensional_gibbs_free_energy_per_link(&1.0, &1.0, nondimensional_link_stiffness, &ZERO, &300.0)
}

/// The implemented functionality of the log-squared-FJC model thermodynamics in the isotensional ensemble.
impl LOGSQUAREDFJC
{
    /// Initializes and returns an instance of the log-squared-FJC model thermodynamics in the isotensional ensemble.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        LOGSQUAREDFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            asymptotic: self::asymptotic::LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness),
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
