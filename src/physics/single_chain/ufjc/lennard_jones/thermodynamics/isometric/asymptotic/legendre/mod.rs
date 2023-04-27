#[cfg(feature = "extern")]
pub mod ex;

#[cfg(feature = "python")]
pub mod py;

mod test;

use std::f64::consts::PI;
use crate::math::
{
    inverse_langevin,
    inverse_newton_raphson_powered
};
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT,
    single_chain::ZERO
};
use crate::physics::single_chain::ufjc::lennard_jones::thermodynamics::nondimensional_link_stretch;

/// The structure of the Lennard-Jones-FJC model thermodynamics in the isometric ensemble approximated using an asymptotic approach and a Legendre transformation.
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

/// The expected force as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, and link stiffness.
pub fn force(number_of_links: &u8, link_length: &f64, link_stiffness: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force(&(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(end_to_end_length/((*number_of_links as f64)*link_length)))
}

/// The expected nondimensional force as a function of the applied nondimensional end-to-end length per link, parameterized by the nondimensional link stiffness.
pub fn nondimensional_force(nondimensional_link_stiffness: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    let c: f64 = 2.0/23.0;
    let lambda_max = (13.0/7.0_f64).powf(1.0/6.0);
    let nondimensional_force_max = nondimensional_link_stiffness/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
    let mut guess: f64 = if nondimensional_end_to_end_length_per_link < &1.0
    {
        inverse_langevin(nondimensional_end_to_end_length_per_link)
    }
    else
    {
        0.95*nondimensional_force_max
    };
    if guess > nondimensional_force_max
    {
        guess = 0.95*nondimensional_force_max;
    }
    inverse_newton_raphson_powered(nondimensional_end_to_end_length_per_link, &|nondimensional_force: &f64| 1.0/nondimensional_force.tanh() - 1.0/nondimensional_force + nondimensional_force/nondimensional_link_stiffness*((nondimensional_force.tanh() - 1.0/nondimensional_force.tanh() + 1.0/nondimensional_force)/(c*nondimensional_force.tanh() + nondimensional_force/nondimensional_link_stiffness)) + nondimensional_link_stretch(nondimensional_link_stiffness, nondimensional_force) - 1.0, &|nondimensional_force: &f64| 1.0/nondimensional_force.powi(2) - 1.0/nondimensional_force.sinh().powi(2) + ((2.0*c*nondimensional_link_stiffness*nondimensional_force/nondimensional_force.tanh() - 2.0*c*nondimensional_link_stiffness + 2.0*nondimensional_force.powi(2) - 1.0)/nondimensional_force.sinh().powi(2) + nondimensional_force.powi(2)/nondimensional_force.sinh().powi(4) - 1.0)/(c*nondimensional_link_stiffness + nondimensional_force/nondimensional_force.tanh()).powi(2) + (nondimensional_link_stretch(nondimensional_link_stiffness, &(nondimensional_force + 1e-6)) - nondimensional_link_stretch(nondimensional_link_stiffness, &(nondimensional_force - 1e-6)))/2e-6, &guess, &1e-6, &100, 4)
}

/// The Helmholtz free energy as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, hinge mass, and link stiffness.
pub fn helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, link_stiffness: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(end_to_end_length/(*number_of_links as f64)/link_length), temperature)
}

/// The Helmholtz free energy per link as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, hinge mass, and link stiffness.
pub fn helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, link_stiffness: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*nondimensional_helmholtz_free_energy_per_link(number_of_links, link_length, hinge_mass, &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(end_to_end_length/(*number_of_links as f64)/link_length), temperature)
}

/// The relative Helmholtz free energy as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, and link stiffness.
pub fn relative_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, link_stiffness: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    helmholtz_free_energy(number_of_links, link_length, &1.0, link_stiffness, end_to_end_length, temperature) - helmholtz_free_energy(number_of_links, link_length, &1.0, link_stiffness, &(ZERO*(*number_of_links as f64)*link_length), temperature)
}

/// The relative Helmholtz free energy per link as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, and link stiffness.
pub fn relative_helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, link_stiffness: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    helmholtz_free_energy_per_link(number_of_links, link_length, &1.0, link_stiffness, end_to_end_length, temperature) - helmholtz_free_energy_per_link(number_of_links, link_length, &1.0, link_stiffness, &(ZERO*(*number_of_links as f64)*link_length), temperature)
}

/// The nondimensional Helmholtz free energy as a function of the applied nondimensional end-to-end length per link and temperature, parameterized by the number of links, link length, hinge mass, and nondimensional link stiffness.
pub fn nondimensional_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_link_stiffness: &f64, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
{
    (*number_of_links as f64)*nondimensional_helmholtz_free_energy_per_link(number_of_links, link_length, hinge_mass, nondimensional_link_stiffness, nondimensional_end_to_end_length_per_link, temperature)
}

/// The nondimensional Helmholtz free energy per link as a function of the nondimensional end-to-end length per link and temperature, parameterized by the number of links, link length, hinge mass, and nondimensional link stiffness.
pub fn nondimensional_helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_link_stiffness: &f64, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
{
    let nondimensional_force = nondimensional_force(nondimensional_link_stiffness, nondimensional_end_to_end_length_per_link);
    let lambda = nondimensional_link_stretch(nondimensional_link_stiffness, &nondimensional_force);
    -(nondimensional_force.sinh()/nondimensional_force).ln() - (1.0 + 23.0/2.0*nondimensional_force/nondimensional_force.tanh()/nondimensional_link_stiffness).ln() + nondimensional_link_stiffness/72.0*(lambda.powi(-12) - 2.0*lambda.powi(-6)) - nondimensional_force*(lambda - 1.0) + nondimensional_force*nondimensional_end_to_end_length_per_link - (1.0 - 1.0/(*number_of_links as f64))*(0.5*(2.0*PI*link_length.powi(2)/nondimensional_link_stiffness).ln() + (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())
}

/// The nondimensional relative Helmholtz free energy as a function of the nondimensional end-to-end length per link, parameterized by the number of links and nondimensional link stiffness.
pub fn nondimensional_relative_helmholtz_free_energy(number_of_links: &u8, nondimensional_link_stiffness: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, &1.0, &1.0, nondimensional_link_stiffness, nondimensional_end_to_end_length_per_link, &300.0) - nondimensional_helmholtz_free_energy(number_of_links, &1.0, &1.0, nondimensional_link_stiffness, &ZERO, &300.0)
}

/// The nondimensional relative Helmholtz free energy per link as a function of the nondimensional end-to-end length per link, parameterized by the nondimensional link stiffness.
pub fn nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_link_stiffness: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    nondimensional_helmholtz_free_energy_per_link(&8, &1.0, &1.0, nondimensional_link_stiffness, nondimensional_end_to_end_length_per_link, &300.0) - nondimensional_helmholtz_free_energy_per_link(&8, &1.0, &1.0, nondimensional_link_stiffness, &ZERO, &300.0)
}

/// The implemented functionality of the Lennard-Jones-FJC model thermodynamics in the isometric ensemble approximated using an asymptotic approach and a Legendre transformation.
impl LENNARDJONESFJC
{
    /// Initializes and returns an instance of the Lennard-Jones-FJC model thermodynamics in the isometric ensemble approximated using an asymptotic approach and a Legendre transformation.
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
    /// The expected force as a function of the applied end-to-end length and temperature.
    pub fn force(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        force(&self.number_of_links, &self.link_length, &self.link_stiffness, end_to_end_length, temperature)
    }
    /// The expected nondimensional force as a function of the applied nondimensional end-to-end length per link and temperature.
    pub fn nondimensional_force(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        nondimensional_force(&(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), nondimensional_end_to_end_length_per_link)
    }
    /// The Helmholtz free energy as a function of the applied end-to-end length and temperature.
    pub fn helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.link_stiffness, end_to_end_length, temperature)
    }
    /// The Helmholtz free energy per link as a function of the applied end-to-end length and temperature.
    pub fn helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.link_stiffness, end_to_end_length, temperature)
    }
    /// The relative Helmholtz free energy as a function of the applied end-to-end length and temperature.
    pub fn relative_helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.link_stiffness, end_to_end_length, temperature)
    }
    /// The relative Helmholtz free energy per link as a function of the applied end-to-end length and temperature.
    pub fn relative_helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.link_stiffness, end_to_end_length, temperature)
    }
    /// The nondimensional Helmholtz free energy as a function of the applied nondimensional end-to-end length per link and temperature.
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), nondimensional_end_to_end_length_per_link, temperature)
    }
    /// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional end-to-end length per link and temperature.
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), nondimensional_end_to_end_length_per_link, temperature)
    }
    /// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy(&self.number_of_links, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), nondimensional_end_to_end_length_per_link)
    }
    /// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy_per_link(&(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), nondimensional_end_to_end_length_per_link)
    }
}
