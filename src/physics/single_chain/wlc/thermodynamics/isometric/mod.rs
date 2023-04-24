#[cfg(feature = "extern")]
pub mod ex;

#[cfg(feature = "python")]
pub mod py;

mod test;

/// The worm-like chain (WLC) model thermodynamics in the isometric ensemble approximated using a Legendre transformation.
pub mod legendre;

use std::f64::consts::PI;
use crate::math::bessel_i;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::ZERO;

/// The structure of the thermodynamics of the WLC model in the isometric ensemble.
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

    contour_length: f64,

    nondimensional_persistance_length: f64,

    /// The thermodynamic functions of the model in the isometric ensemble approximated using a Legendre transformation.
    pub legendre: self::legendre::WLC
}

/// The expected force as a function of the applied end-to-end length and temperature, parameterized by the link length, contour length, and persistance length.
pub fn force(link_length: &f64, contour_length: &f64, persistance_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    nondimensional_force(&(persistance_length/contour_length), &(end_to_end_length/contour_length))*BOLTZMANN_CONSTANT*temperature/link_length
}

/// The expected nondimensional force as a function of the applied nondimensional end-to-end length per link, parameterized by the nondimensional persistance length.
pub fn nondimensional_force(nondimensional_persistance_length: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    bessel_i(&1, nondimensional_end_to_end_length_per_link)*nondimensional_persistance_length
}

//
//
//
// use logarithm of distributions for free energies?
//
//
//

/// The equilibrium probability density of end-to-end vectors as a function of the end-to-end length, parameterized by the contour length and persistance length.
pub fn equilibrium_distribution(contour_length: &f64, persistance_length: &f64, end_to_end_length: &f64) -> f64
{
    nondimensional_equilibrium_distribution(&(persistance_length/contour_length), &(end_to_end_length/contour_length))/contour_length.powi(3)
}

/// The nondimensional equilibrium probability density of nondimensional end-to-end vectors per link as a function of the nondimensional end-to-end length per link, parameterized by the nondimensional persistance length.
pub fn nondimensional_equilibrium_distribution(nondimensional_persistance_length: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    bessel_i(&0, nondimensional_end_to_end_length_per_link)*nondimensional_persistance_length
}

/// The equilibrium probability density of end-to-end lengths as a function of the end-to-end length, parameterized by the contour length and persistance length.
pub fn equilibrium_radial_distribution(contour_length: &f64, persistance_length: &f64, end_to_end_length: &f64) -> f64
{
    nondimensional_equilibrium_radial_distribution(&(persistance_length/contour_length), &(end_to_end_length/contour_length))/contour_length
}

/// The nondimensional equilibrium probability density of nondimensional end-to-end lengths per link as a function of the nondimensional end-to-end length per link, parameterized by the nondimensional persistance length.
pub fn nondimensional_equilibrium_radial_distribution(nondimensional_persistance_length: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    4.0*PI*nondimensional_end_to_end_length_per_link.powi(2)*nondimensional_equilibrium_distribution(nondimensional_persistance_length, nondimensional_end_to_end_length_per_link)
}

/// The implemented functionality of the thermodynamics of the WLC model in the isometric ensemble.
impl WLC
{
    /// Initializes and returns an instance of the thermodynamics of the WLC model in the isometric ensemble.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, persistance_length: f64) -> Self
    {
        WLC
        {
            hinge_mass,
            link_length,
            number_of_links,
            persistance_length,
            contour_length: (number_of_links as f64)*link_length,
            nondimensional_persistance_length: persistance_length/(number_of_links as f64)/link_length,
            legendre: self::legendre::WLC::init(number_of_links, link_length, hinge_mass, persistance_length),
        }
    }
    /// The expected force as a function of the applied end-to-end length and temperature.
    pub fn force(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        force(&self.link_length, &self.contour_length, &self.persistance_length, end_to_end_length, temperature)
    }
    /// The expected nondimensional force as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_force(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_force(&self.nondimensional_persistance_length, nondimensional_end_to_end_length_per_link)
    }
    /// The equilibrium probability density of end-to-end vectors as a function of the end-to-end length.
    pub fn equilibrium_distribution(&self, end_to_end_length: &f64) -> f64
    {
        equilibrium_distribution(&self.contour_length, &self.persistance_length, end_to_end_length)
    }
    /// The nondimensional equilibrium probability density of nondimensional end-to-end vectors per link as a function of the nondimensional end-to-end length per link.
    pub fn nondimensional_equilibrium_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_equilibrium_distribution(&self.nondimensional_persistance_length, nondimensional_end_to_end_length_per_link)
    }
    /// The equilibrium probability density of end-to-end lengths as a function of the end-to-end length.
    pub fn equilibrium_radial_distribution(&self, end_to_end_length: &f64) -> f64
    {
        equilibrium_radial_distribution(&self.contour_length, &self.persistance_length, end_to_end_length)
    }
    /// The nondimensional equilibrium probability density of nondimensional end-to-end lengths per link as a function of the nondimensional end-to-end length per link.
    pub fn nondimensional_equilibrium_radial_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_equilibrium_radial_distribution(&self.nondimensional_persistance_length, nondimensional_end_to_end_length_per_link)
    }
}
