#[cfg(feature = "extern")]
pub mod ex;

#[cfg(feature = "python")]
pub mod py;

mod test;

/// The freely-jointed chain (FJC) model thermodynamics in the isometric ensemble approximated using a Legendre transformation.
pub mod legendre;

use super::treloar_sums;
use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::ZERO;

/// The structure of the thermodynamics of the FJC model in the isometric ensemble.
pub struct FJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The thermodynamic functions of the model in the isometric ensemble approximated using a Legendre transformation.
    pub legendre: legendre::FJC
}

/// The expected force as a function of the applied end-to-end length and temperature, parameterized by the number of links and link length.
pub fn force(number_of_links: &u8, link_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    nondimensional_force(number_of_links, &(end_to_end_length/((*number_of_links as f64)*link_length)))*BOLTZMANN_CONSTANT*temperature/link_length
}

/// The expected nondimensional force as a function of the applied nondimensional end-to-end length per link, parameterized by the number of links.
pub fn nondimensional_force(number_of_links: &u8, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    let sums = treloar_sums(number_of_links, nondimensional_end_to_end_length_per_link, &vec![0, 1]);
    let number_of_links_f64 = *number_of_links as f64;
    (1.0/nondimensional_end_to_end_length_per_link + (0.5*number_of_links_f64 - 1.0)*sums[1]/sums[0])/number_of_links_f64
}

/// The Helmholtz free energy as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, &(end_to_end_length/((*number_of_links as f64)*link_length)), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The Helmholtz free energy per link as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy_per_link(number_of_links, link_length, hinge_mass, &(end_to_end_length/((*number_of_links as f64)*link_length)), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The relative Helmholtz free energy as a function of the applied end-to-end length and temperature, parameterized by the number of links and link length.
pub fn relative_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy(number_of_links, &(end_to_end_length/((*number_of_links as f64)*link_length)))*BOLTZMANN_CONSTANT*temperature
}

/// The relative Helmholtz free energy per link as a function of the applied end-to-end length and temperature, parameterized by the number of links and link length.
pub fn relative_helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy_per_link(number_of_links, &(end_to_end_length/((*number_of_links as f64)*link_length)))*BOLTZMANN_CONSTANT*temperature
}

/// The nondimensional Helmholtz free energy as a function of the nondimensional end-to-end length per link and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
{
    -(equilibrium_distribution(number_of_links, link_length, &(nondimensional_end_to_end_length_per_link*(*number_of_links as f64)*link_length))).ln() - ((*number_of_links as f64) - 1.0)*(8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional Helmholtz free energy per link as a function of the nondimensional end-to-end length per link and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, nondimensional_end_to_end_length_per_link, temperature)/(*number_of_links as f64)
}

/// The nondimensional relative Helmholtz free energy as a function of the nondimensional end-to-end length per link, parameterized by the number of links.
pub fn nondimensional_relative_helmholtz_free_energy(number_of_links: &u8, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    (nondimensional_equilibrium_distribution(number_of_links, &ZERO)/nondimensional_equilibrium_distribution(number_of_links, nondimensional_end_to_end_length_per_link)).ln()
}

/// The nondimensional relative Helmholtz free energy per link as a function of the nondimensional end-to-end length per link, parameterized by the number of links.
pub fn nondimensional_relative_helmholtz_free_energy_per_link(number_of_links: &u8, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy(number_of_links, nondimensional_end_to_end_length_per_link)/(*number_of_links as f64)
}

/// The equilibrium probability density of end-to-end vectors as a function of the end-to-end length, parameterized by the number of links and link length.
pub fn equilibrium_distribution(number_of_links: &u8, link_length: &f64, end_to_end_length: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_equilibrium_distribution(number_of_links, &(end_to_end_length/contour_length))/contour_length.powi(3)
}

/// The nondimensional equilibrium probability density of nondimensional end-to-end vectors per link as a function of the nondimensional end-to-end length per link, parameterized by the number of links.
pub fn nondimensional_equilibrium_distribution(number_of_links: &u8, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    let number_of_links_f64 = *number_of_links as f64;
    let n = *number_of_links as u128;
    let p: i32 = (number_of_links - 2).into();
    let m = -*nondimensional_end_to_end_length_per_link*0.5 + 0.5;
    let k = (number_of_links_f64*m).ceil() as u128;
    let sum: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/number_of_links_f64).powi(p)).sum();
    0.125/PI/nondimensional_end_to_end_length_per_link*(n.pow(n as u32) as f64)/((1..=n-2).product::<u128>() as f64)*sum
}

/// The equilibrium probability density of end-to-end lengths as a function of the end-to-end length, parameterized by the number of links and link length.
pub fn equilibrium_radial_distribution(number_of_links: &u8, link_length: &f64, end_to_end_length: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_equilibrium_radial_distribution(number_of_links, &(end_to_end_length/contour_length))/contour_length
}

/// The nondimensional equilibrium probability density of nondimensional end-to-end lengths per link as a function of the nondimensional end-to-end length per link, parameterized by the number of links.
pub fn nondimensional_equilibrium_radial_distribution(number_of_links: &u8, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    let number_of_links_f64 = *number_of_links as f64;
    let n = *number_of_links as u128;
    let p: i32 = (number_of_links - 2).into();
    let m = -*nondimensional_end_to_end_length_per_link*0.5 + 0.5;
    let k = (number_of_links_f64*m).ceil() as u128;
    let sum: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/number_of_links_f64).powi(p)).sum();
    0.5*nondimensional_end_to_end_length_per_link*(n.pow(n as u32) as f64)/((1..=n-2).product::<u128>() as f64)*sum
}

/// The implemented functionality of the thermodynamics of the FJC model in the isometric ensemble.
impl FJC
{
    /// Initializes and returns an instance of the thermodynamics of the FJC model in the isometric ensemble.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            legendre: legendre::FJC::init(number_of_links, link_length, hinge_mass)
        }
    }
    /// The expected force as a function of the applied end-to-end length and temperature.
    pub fn force(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        force(&self.number_of_links, &self.link_length, end_to_end_length, temperature)
    }
    /// The expected nondimensional force as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_force(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_force(&self.number_of_links, nondimensional_end_to_end_length_per_link)
    }
    /// The Helmholtz free energy as a function of the applied end-to-end length and temperature.
    pub fn helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, end_to_end_length, temperature)
    }
    /// The Helmholtz free energy per link as a function of the applied end-to-end length and temperature.
    pub fn helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, end_to_end_length, temperature)
    }
    /// The relative Helmholtz free energy as a function of the applied end-to-end length and temperature.
    pub fn relative_helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy(&self.number_of_links, &self.link_length, end_to_end_length, temperature)
    }
    /// The relative Helmholtz free energy per link as a function of the applied end-to-end length and temperature.
    pub fn relative_helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, end_to_end_length, temperature)
    }
    /// The nondimensional Helmholtz free energy as a function of the applied nondimensional end-to-end length per link and temperature.
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, nondimensional_end_to_end_length_per_link, temperature)
    }
    /// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional end-to-end length per link and temperature.
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, nondimensional_end_to_end_length_per_link, temperature)
    }
    /// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy(&self.number_of_links, nondimensional_end_to_end_length_per_link)
    }
    /// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy_per_link(&self.number_of_links, nondimensional_end_to_end_length_per_link)
    }
    /// The equilibrium probability density of end-to-end vectors as a function of the end-to-end length.
    pub fn equilibrium_distribution(&self, end_to_end_length: &f64) -> f64
    {
        equilibrium_distribution(&self.number_of_links, &self.link_length, end_to_end_length)
    }
    /// The nondimensional equilibrium probability density of nondimensional end-to-end vectors per link as a function of the nondimensional end-to-end length per link.
    pub fn nondimensional_equilibrium_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_equilibrium_distribution(&self.number_of_links, nondimensional_end_to_end_length_per_link)
    }
    /// The equilibrium probability density of end-to-end lengths as a function of the end-to-end length.
    pub fn equilibrium_radial_distribution(&self, end_to_end_length: &f64) -> f64
    {
        equilibrium_radial_distribution(&self.number_of_links, &self.link_length, end_to_end_length)
    }
    /// The nondimensional equilibrium probability density of nondimensional end-to-end lengths per link as a function of the nondimensional end-to-end length per link.
    pub fn nondimensional_equilibrium_radial_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_equilibrium_radial_distribution(&self.number_of_links, nondimensional_end_to_end_length_per_link)
    }
}
