#[cfg(feature = "python")]
pub mod py;

mod test;

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

/// The structure of the Lennard-Jones-FJC model thermodynamics in the isotensional ensemble approximated using a Legendre transformation.
pub struct LENNARDJONESFJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The stiffness of each link in the chain in units of J/(mol⋅nm^2).
    pub link_stiffness: f64
}

/// The implemented functionality of the Lennard-Jones-FJC model thermodynamics in the isotensional ensemble approximated using a Legendre transformation.
impl LENNARDJONESFJC
{
    /// Initializes and returns an instance of the Lennard-Jones-FJC model thermodynamics in the isotensional ensemble approximated using a Legendre transformation.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        LENNARDJONESFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness
        }
    }
}
