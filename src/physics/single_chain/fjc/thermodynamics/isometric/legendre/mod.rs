#[cfg(feature = "extern")]
pub mod ex;

#[cfg(feature = "python")]
pub mod py;

mod test;

use super::
{
    treloar_sums,
    treloar_sum_0_with_prefactor
};
use std::f64::consts::PI;
use crate::math::
{
    inverse_langevin,
    integrate_1d
};
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

/// The structure of the thermodynamics of the FJC model in the isometric ensemble approximated using a Legendre transformation.
pub struct FJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,
    
    normalization_nondimensional_equilibrium_distribution: f64
}

/// The expected force as a function of the applied end-to-end length and temperature, parameterized by the number of links and link length.
pub fn force(number_of_links: &u8, link_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force(&(end_to_end_length/((*number_of_links as f64)*link_length)))
}

/// The expected nondimensional force as a function of the applied nondimensional end-to-end length per link, parameterized by the number of links.
pub fn nondimensional_force(nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    inverse_langevin(nondimensional_end_to_end_length_per_link)
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
    nondimensional_relative_helmholtz_free_energy_per_link(&(end_to_end_length/((*number_of_links as f64)*link_length)))*BOLTZMANN_CONSTANT*temperature
}

/// The nondimensional Helmholtz free energy as a function of the nondimensional end-to-end length per link and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy(number_of_links, nondimensional_end_to_end_length_per_link) - ((*number_of_links as f64) - 1.0)*(8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional Helmholtz free energy per link as a function of the nondimensional end-to-end length per link and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_end_to_end_length_per_link) - (1.0 - 1.0/(*number_of_links as f64))*(8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional relative Helmholtz free energy as a function of the nondimensional end-to-end length per link, parameterized by the number of links.
pub fn nondimensional_relative_helmholtz_free_energy(number_of_links: &u8, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_end_to_end_length_per_link)*(*number_of_links as f64)
}

/// The nondimensional relative Helmholtz free energy per link as a function of the nondimensional end-to-end length per link, parameterized by the number of links.
pub fn nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    let nondimensional_force = nondimensional_force(nondimensional_end_to_end_length_per_link);
    nondimensional_force**nondimensional_end_to_end_length_per_link - (nondimensional_force.sinh()/nondimensional_force).ln()
}

/// The equilibrium probability density of end-to-end vectors as a function of the end-to-end length, parameterized by the number of links and link length.
pub fn equilibrium_distribution(number_of_links: &u8, link_length: &f64, normalization_nondimensional_equilibrium_distribution: &f64, end_to_end_length: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_equilibrium_distribution(number_of_links, normalization_nondimensional_equilibrium_distribution, &(end_to_end_length/contour_length))/contour_length.powi(3)
}

/// The nondimensional equilibrium probability density of nondimensional end-to-end vectors per link as a function of the nondimensional end-to-end length per link, parameterized by the number of links.
pub fn nondimensional_equilibrium_distribution(number_of_links: &u8, normalization_nondimensional_equilibrium_distribution: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    let nondimensional_force = nondimensional_force(nondimensional_end_to_end_length_per_link);
(nondimensional_force.sinh()/nondimensional_force*(-nondimensional_force**nondimensional_end_to_end_length_per_link).exp()).powi(*number_of_links as i32)/normalization_nondimensional_equilibrium_distribution
}

/// The equilibrium probability density of end-to-end lengths as a function of the end-to-end length, parameterized by the number of links and link length.
pub fn equilibrium_radial_distribution(number_of_links: &u8, link_length: &f64, normalization_nondimensional_equilibrium_distribution: &f64, end_to_end_length: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_equilibrium_radial_distribution(number_of_links, normalization_nondimensional_equilibrium_distribution, &(end_to_end_length/contour_length))/contour_length
}

/// The nondimensional equilibrium probability density of nondimensional end-to-end lengths per link as a function of the nondimensional end-to-end length per link, parameterized by the number of links.
pub fn nondimensional_equilibrium_radial_distribution(number_of_links: &u8, normalization_nondimensional_equilibrium_distribution: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    4.0*PI*nondimensional_end_to_end_length_per_link.powi(2)*nondimensional_equilibrium_distribution(number_of_links, normalization_nondimensional_equilibrium_distribution, nondimensional_end_to_end_length_per_link)
}

/// The Gibbs free energy as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    nondimensional_gibbs_free_energy(number_of_links, link_length, hinge_mass, &(end_to_end_length/((*number_of_links as f64)*link_length)), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The Gibbs free energy per link as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    nondimensional_gibbs_free_energy_per_link(number_of_links, link_length, hinge_mass, &(end_to_end_length/((*number_of_links as f64)*link_length)), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The relative Gibbs free energy as a function of the applied end-to-end length and temperature, parameterized by the number of links and link length.
pub fn relative_gibbs_free_energy(number_of_links: &u8, link_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    nondimensional_relative_gibbs_free_energy(number_of_links, &(end_to_end_length/((*number_of_links as f64)*link_length)))*BOLTZMANN_CONSTANT*temperature
}

/// The relative Gibbs free energy per link as a function of the applied end-to-end length and temperature, parameterized by the number of links and link length.
pub fn relative_gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    nondimensional_relative_gibbs_free_energy_per_link(number_of_links, &(end_to_end_length/((*number_of_links as f64)*link_length)))*BOLTZMANN_CONSTANT*temperature
}

/// The nondimensional Gibbs free energy as a function of the applied nondimensional end-to-end length per link and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
{
    let sums = treloar_sums(number_of_links, nondimensional_end_to_end_length_per_link, &[0, 1]);
    let number_of_links_f64 = *number_of_links as f64;
    -(nondimensional_end_to_end_length_per_link*number_of_links_f64)*(1.0/nondimensional_end_to_end_length_per_link + (0.5*number_of_links_f64 - 1.0)*sums[1]/sums[0])/number_of_links_f64 - (treloar_sum_0_with_prefactor(number_of_links, nondimensional_end_to_end_length_per_link)/(number_of_links_f64*link_length).powi(3)).ln() - (number_of_links_f64 - 1.0)*(8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional Gibbs free energy per link as a function of the applied nondimensional end-to-end length per link and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn nondimensional_gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
{
    nondimensional_gibbs_free_energy(number_of_links, link_length, hinge_mass, nondimensional_end_to_end_length_per_link, temperature)/(*number_of_links as f64)
}

/// The nondimensional relative Gibbs free energy as a function of the applied nondimensional end-to-end length per link, parameterized by the number of links.
pub fn nondimensional_relative_gibbs_free_energy(number_of_links: &u8, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    nondimensional_gibbs_free_energy(number_of_links, &1.0, &1.0, nondimensional_end_to_end_length_per_link, &300.0) - nondimensional_gibbs_free_energy(number_of_links, &1.0, &1.0, &ZERO, &300.0)
}

/// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional end-to-end length per link, parameterized by the number of links.
pub fn nondimensional_relative_gibbs_free_energy_per_link(number_of_links: &u8, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    nondimensional_gibbs_free_energy_per_link(number_of_links, &1.0, &1.0, nondimensional_end_to_end_length_per_link, &300.0) - nondimensional_gibbs_free_energy_per_link(number_of_links, &1.0, &1.0, &ZERO, &300.0)
}

/// The implemented functionality of the thermodynamics of the FJC model in the isometric ensemble approximated using a Legendre transformation.
impl FJC
{
    /// Initializes and returns an instance of the thermodynamics of the FJC model in the isometric ensemble approximated using a Legendre transformation.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        let normalization = integrate_1d(&|nondimensional_end_to_end_length_per_link| nondimensional_equilibrium_radial_distribution(&number_of_links, &1.0, nondimensional_end_to_end_length_per_link), &ZERO, &ONE, &POINTS);
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            normalization_nondimensional_equilibrium_distribution: normalization
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
        nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_end_to_end_length_per_link)
    }
    /// The equilibrium probability density of end-to-end vectors as a function of the end-to-end length.
    pub fn equilibrium_distribution(&self, end_to_end_length: &f64) -> f64
    {
        equilibrium_distribution(&self.number_of_links, &self.link_length, &self.normalization_nondimensional_equilibrium_distribution, end_to_end_length)
    }
    /// The equilibrium probability density of nondimensional end-to-end vectors per link as a function of the nondimensional end-to-end length per link.
    pub fn nondimensional_equilibrium_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_equilibrium_distribution(&self.number_of_links, &self.normalization_nondimensional_equilibrium_distribution, nondimensional_end_to_end_length_per_link)
    }
    /// The equilibrium probability density of end-to-end lengths as a function of the end-to-end length.
    pub fn equilibrium_radial_distribution(&self, end_to_end_length: &f64) -> f64
    {
        equilibrium_radial_distribution(&self.number_of_links, &self.link_length, &self.normalization_nondimensional_equilibrium_distribution, end_to_end_length)
    }
    /// The equilibrium probability density of nondimensional end-to-end lengths per link as a function of the nondimensional end-to-end length per link.
    pub fn nondimensional_equilibrium_radial_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_equilibrium_radial_distribution(&self.number_of_links, &self.normalization_nondimensional_equilibrium_distribution, nondimensional_end_to_end_length_per_link)
    }
    /// The Gibbs free energy as a function of the applied end-to-end length and temperature.
    pub fn gibbs_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, end_to_end_length, temperature)
    }
    /// The Gibbs free energy per link as a function of the applied end-to-end length and temperature.
    pub fn gibbs_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, end_to_end_length, temperature)
    }
    /// The relative Gibbs free energy as a function of the applied end-to-end length and temperature.
    pub fn relative_gibbs_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        relative_gibbs_free_energy(&self.number_of_links, &self.link_length, end_to_end_length, temperature)
    }
    /// The relative Gibbs free energy per link as a function of the applied end-to-end length and temperature.
    pub fn relative_gibbs_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        relative_gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, end_to_end_length, temperature)
    }
    /// The nondimensional Gibbs free energy as a function of the applied nondimensional end-to-end length per link and temperature.
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        nondimensional_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, nondimensional_end_to_end_length_per_link, temperature)
    }
    /// The nondimensional Gibbs free energy per link as a function of the applied nondimensional end-to-end length per link and temperature.
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        nondimensional_gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, nondimensional_end_to_end_length_per_link, temperature)
    }
    /// The nondimensional relative Gibbs free energy as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_relative_gibbs_free_energy(&self.number_of_links, nondimensional_end_to_end_length_per_link)
    }
    /// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_relative_gibbs_free_energy_per_link(&self.number_of_links, nondimensional_end_to_end_length_per_link)
    }
}
