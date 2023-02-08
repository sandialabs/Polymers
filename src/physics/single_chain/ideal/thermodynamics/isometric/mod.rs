#[cfg(feature = "extern")]
pub mod ex;

#[cfg(feature = "python")]
pub mod py;

mod test;

use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};

/// The structure of the thermodynamics of the ideal chain model in the isometric ensemble.
pub struct Ideal
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8
}

fn force(number_of_links: &u8, link_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    3.0*end_to_end_length*BOLTZMANN_CONSTANT*temperature/(*number_of_links as f64)/link_length.powi(2)
}

fn nondimensional_force(nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    3.0*nondimensional_end_to_end_length_per_link
}

fn helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    (*number_of_links as f64)*BOLTZMANN_CONSTANT*temperature*(1.5*(end_to_end_length/((*number_of_links as f64)*link_length)).powi(2) - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())
}

fn helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*(1.5*(end_to_end_length/((*number_of_links as f64)*link_length)).powi(2) - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())
}

fn relative_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy(number_of_links, &(end_to_end_length/((*number_of_links as f64)*link_length)))*BOLTZMANN_CONSTANT*temperature
}

fn relative_helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy_per_link(&(end_to_end_length/((*number_of_links as f64)*link_length)))*BOLTZMANN_CONSTANT*temperature
}

fn nondimensional_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
{
    (*number_of_links as f64)*(1.5*nondimensional_end_to_end_length_per_link.powi(2) - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())
}

fn nondimensional_helmholtz_free_energy_per_link(link_length: &f64, hinge_mass: &f64, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
{
    1.5*nondimensional_end_to_end_length_per_link.powi(2) - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

fn nondimensional_relative_helmholtz_free_energy(number_of_links: &u8, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    1.5*(*number_of_links as f64)*nondimensional_end_to_end_length_per_link.powi(2)
}

fn nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    1.5*nondimensional_end_to_end_length_per_link.powi(2)
}

fn equilibrium_distribution(number_of_links: &u8, link_length: &f64, end_to_end_length: &f64) -> f64
{
    (1.5/PI/(*number_of_links as f64)/link_length.powi(2)).powf(1.5)*(-1.5*(end_to_end_length/link_length).powi(2)/(*number_of_links as f64)).exp()
}

fn nondimensional_equilibrium_distribution(number_of_links: &u8, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    (1.5/PI*(*number_of_links as f64)).powf(1.5)*(-1.5*nondimensional_end_to_end_length_per_link.powi(2)*(*number_of_links as f64)).exp()
}

fn equilibrium_radial_distribution(number_of_links: &u8, link_length: &f64, end_to_end_length: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    4.0*PI*(end_to_end_length/contour_length).powi(2)*(1.5/PI*(*number_of_links as f64)).powf(1.5)*(-1.5*(end_to_end_length/contour_length).powi(2)*(*number_of_links as f64)).exp()/contour_length
}

fn nondimensional_equilibrium_radial_distribution(number_of_links: &u8, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    4.0*PI*nondimensional_end_to_end_length_per_link.powi(2)*(1.5/PI*(*number_of_links as f64)).powf(1.5)*(-1.5*nondimensional_end_to_end_length_per_link.powi(2)*(*number_of_links as f64)).exp()
}

/// The implemented functionality of the thermodynamics of the ideal chain model in the isometric ensemble.
impl Ideal
{
    /// Initializes and returns an instance of the thermodynamics of the ideal chain model in the isometric ensemble.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        Ideal
        {
            hinge_mass,
            link_length,
            number_of_links
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
        nondimensional_force(nondimensional_end_to_end_length_per_link)
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
        nondimensional_helmholtz_free_energy_per_link(&self.link_length, &self.hinge_mass, nondimensional_end_to_end_length_per_link, temperature)
    }
    /// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy(&self.number_of_links, nondimensional_end_to_end_length_per_link)
    }
    /// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_end_to_end_length_per_link)
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
