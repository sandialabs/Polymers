#[cfg(feature = "python")]
pub mod py;

mod test;

use crate::math::inverse_newton_raphson;

/// The Lennard-Jones link potential freely-jointed chain (Lennard-Jones-FJC) model thermodynamics in the isometric ensemble.
pub mod isometric;

/// The Lennard-Jones link potential freely-jointed chain (Lennard-Jones-FJC) model thermodynamics in the isotensional ensemble.
pub mod isotensional;

/// The structure of the Lennard-Jones-FJC model thermodynamics.
pub struct LENNARDJONESFJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The stiffness of each link in the chain in units of J/(mol⋅nm^2).
    pub link_stiffness: f64,

    /// The thermodynamic functions of the model in the isometric ensemble.
    pub isometric: self::isometric::LENNARDJONESFJC,

    /// The thermodynamic functions of the model in the isotensional ensemble.
    pub isotensional: self::isotensional::LENNARDJONESFJC
}

/// The implemented functionality of the Lennard-Jones-FJC model thermodynamics.
impl LENNARDJONESFJC
{
    /// Initializes and returns an instance of the Lennard-Jones-FJC model thermodynamics.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        LENNARDJONESFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            isometric: self::isometric::LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness),
            isotensional: self::isotensional::LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness)
        }
    }
}

fn nondimensional_link_stretch(nondimensional_link_stiffness: &f64, nondimensional_force: &f64) -> f64
{
    inverse_newton_raphson(&(6.0*nondimensional_force/nondimensional_link_stiffness), &|x: &f64| x.powi(-7) - x.powi(-13), &|x: &f64| 13.0*x.powi(-14) - 7.0*x.powi(-8), &1.0, &1e-6, &100)
}
