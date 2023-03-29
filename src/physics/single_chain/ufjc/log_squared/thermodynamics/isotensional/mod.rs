#[cfg(feature = "python")]
pub mod py;

mod test;

/// The log-squared link potential freely-jointed chain (log-squared-FJC) model thermodynamics in the isotensional ensemble approximated using an asymptotic approach.
pub mod asymptotic;

/// The log-squared link potential freely-jointed chain (log-squared-FJC) model thermodynamics in the isotensional ensemble approximated using a Legendre transformation.
pub mod legendre;

use super::lambert_w;
use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT,
    single_chain::
    {
        ONE,
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
}
