#[cfg(feature = "extern")]
pub mod ex;

#[cfg(feature = "python")]
pub mod py;

mod test;

use super::super::super::
{
    treloar_sums,
    treloar_sum_0_with_prefactor
};
use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::ZERO;

/// The structure of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach valid for strong potentials.
pub struct FJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8
}

/// The expected force as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links and link length.
pub fn force(number_of_links: &u8, link_length: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force(number_of_links, &(*potential_distance/((*number_of_links as f64)*link_length)), &(potential_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature))
}

/// The expected nondimensional force as a function of the applied nondimensional potential distance and nondimensional potential stiffness, parameterized by the number of links.
pub fn nondimensional_force(number_of_links: &u8, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
{
    let sums = treloar_sums(number_of_links, nondimensional_potential_distance, &[0, 1, 2, 3]);
    let number_of_links_f64 = *number_of_links as f64;
    (1.0/nondimensional_potential_distance + (0.5*number_of_links_f64 - 1.0)*sums[1]/sums[0])/number_of_links_f64 + 0.5/nondimensional_potential_stiffness/number_of_links_f64.powi(3)*((0.5*number_of_links_f64 - 1.0)*((0.5*number_of_links_f64 - 1.0)*sums[1]/sums[0]*((number_of_links_f64 - 2.0)*(sums[1]/sums[0]).powi(2) - (number_of_links_f64 - 3.0)*sums[2]/sums[0]) - (0.5*number_of_links_f64 - 1.5)*((0.5*number_of_links_f64 - 1.0)*sums[1]*sums[2]/sums[0].powi(2) - (0.5*number_of_links_f64 - 2.0)*sums[3]/sums[0])) + 2.0*nondimensional_potential_distance.powi(-3) - 2.0*((0.5*number_of_links_f64 - 1.0)*sums[1]/sums[0] + nondimensional_potential_distance.powi(-1))*((0.5*number_of_links_f64 - 1.0)*((0.5*number_of_links_f64 - 1.0)*(sums[1]/sums[0]).powi(2) - (0.5*number_of_links_f64 - 1.5)*sums[2]/sums[0]) - nondimensional_potential_distance.powi(-2)))
}

/// The Helmholtz free energy as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, &(potential_distance/((*number_of_links as f64)*link_length)), &(potential_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The Helmholtz free energy per link as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy_per_link(number_of_links, link_length, hinge_mass, &(potential_distance/((*number_of_links as f64)*link_length)), &(potential_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The relative Helmholtz free energy as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links and link length.
pub fn relative_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    helmholtz_free_energy(number_of_links, link_length, &1.0, potential_distance, potential_stiffness, temperature) - helmholtz_free_energy(number_of_links, link_length, &1.0, &(ZERO*(*number_of_links as f64)*link_length), potential_stiffness, temperature)
}

/// The relative Helmholtz free energy per link as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links and link length.
pub fn relative_helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    helmholtz_free_energy_per_link(number_of_links, link_length, &1.0, potential_distance, potential_stiffness, temperature) - helmholtz_free_energy_per_link(number_of_links, link_length, &1.0, &(ZERO*(*number_of_links as f64)*link_length), potential_stiffness, temperature)
}

/// The nondimensional Helmholtz free energy as a function of the applied nondimensional potential distance, nondimensional potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
{
    let sums = treloar_sums(number_of_links, nondimensional_potential_distance, &[0, 1, 2]);
    let number_of_links_f64 = *number_of_links as f64;
    let number_of_links_squared = number_of_links_f64.powi(2);
    let contour_length = number_of_links_f64*link_length;
    -(treloar_sum_0_with_prefactor(number_of_links, nondimensional_potential_distance)/contour_length.powi(3)).ln() - (number_of_links_f64 - 1.0)*(8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln() - 1.5*(2.0*PI/nondimensional_potential_stiffness/number_of_links_squared).ln() - 3.0*(contour_length).ln() + 0.5/nondimensional_potential_stiffness/number_of_links_squared*((0.5*number_of_links_f64 - 1.0)*((0.5*number_of_links_f64 - 1.0)*(sums[1]/sums[0]).powi(2) - (0.5*number_of_links_f64 - 1.5)*sums[2]/sums[0]) - nondimensional_potential_distance.powi(-2) - ((0.5*number_of_links_f64 - 1.0)*sums[1]/sums[0] + nondimensional_potential_distance.powi(-1)).powi(2))
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

/// The implemented functionality of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach valid for strong potentials.
impl FJC
{
    /// Initializes and returns an instance of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach valid for strong potentials.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links
        }
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
}
