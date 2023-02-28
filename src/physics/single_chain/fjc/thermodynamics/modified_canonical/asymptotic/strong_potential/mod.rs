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
    let contour_length = (*number_of_links as f64)*link_length;
    BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force(number_of_links, &(*potential_distance/contour_length), &(potential_stiffness*contour_length.powi(2)/BOLTZMANN_CONSTANT/temperature))
}

/// The expected nondimensional force as a function of the applied nondimensional potential distance and nondimensional potential stiffness, parameterized by the number of links.
pub fn nondimensional_force(number_of_links: &u8, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
{
    let number_of_links_f64 = *number_of_links as f64;
    let n = *number_of_links as u128;
    let p: i32 = (number_of_links - 2).into();
    let m = -*nondimensional_potential_distance*0.5 + 0.5;
    let k = (number_of_links_f64*m).ceil() as u128;
    let sum_0: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/number_of_links_f64).powi(p)).sum();
    let sum_1: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/number_of_links_f64).powi(p - 1)).sum();
    let sum_2: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/number_of_links_f64).powi(p - 2)).sum();
    let sum_3: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/number_of_links_f64).powi(p - 3)).sum();
    (1.0/nondimensional_potential_distance + (0.5*number_of_links_f64 - 1.0)*sum_1/sum_0)/number_of_links_f64 + 0.5/nondimensional_potential_stiffness/number_of_links_f64*((0.5*number_of_links_f64 - 1.0)*((0.5*number_of_links_f64 - 1.0)*sum_1/sum_0*((number_of_links_f64 - 2.0)*(sum_1/sum_0).powi(2) - (number_of_links_f64 - 3.0)*sum_2/sum_0) - (0.5*number_of_links_f64 - 1.5)*((0.5*number_of_links_f64 - 1.0)*sum_1*sum_2/sum_0.powi(2) - (0.5*number_of_links_f64 - 2.0)*sum_3/sum_0)) + 2.0*nondimensional_potential_distance.powi(-3) - 2.0*((0.5*number_of_links_f64 - 1.0)*sum_1/sum_0 + nondimensional_potential_distance.powi(-1))*((0.5*number_of_links_f64 - 1.0)*((0.5*number_of_links_f64 - 1.0)*(sum_1/sum_0).powi(2) - (0.5*number_of_links_f64 - 1.5)*sum_2/sum_0) - nondimensional_potential_distance.powi(-2)))
}

/// The Helmholtz free energy as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, &(potential_distance/contour_length), &(potential_stiffness*contour_length.powi(2)/BOLTZMANN_CONSTANT/temperature), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The Helmholtz free energy per link as a function of the applied potential distance, potential stiffness, and temperature, parameterized by the number of links, link length, and hinge mass.
pub fn helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_helmholtz_free_energy_per_link(number_of_links, link_length, hinge_mass, &(potential_distance/contour_length), &(potential_stiffness*contour_length.powi(2)/BOLTZMANN_CONSTANT/temperature), temperature)*BOLTZMANN_CONSTANT*temperature
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
    let number_of_links_f64 = *number_of_links as f64;
    let contour_length = number_of_links_f64*link_length;
    let n = *number_of_links as u128;
    let p: i32 = (number_of_links - 2).into();
    let m = -*nondimensional_potential_distance*0.5 + 0.5;
    let k = (number_of_links_f64*m).ceil() as u128;
    let sum_0: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/number_of_links_f64).powi(p)).sum();
    let sum_1: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/number_of_links_f64).powi(p - 1)).sum();
    let sum_2: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/number_of_links_f64).powi(p - 2)).sum();
    -(0.125/PI/nondimensional_potential_distance*(n.pow(n as u32) as f64)/((1..=n-2).product::<u128>() as f64)*sum_0/contour_length.powi(3)).ln() - (number_of_links_f64 - 1.0)*(8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln() - 1.5*(2.0*PI/nondimensional_potential_stiffness).ln() - 3.0*(contour_length).ln() + 0.5/nondimensional_potential_stiffness*((0.5*number_of_links_f64 - 1.0)*((0.5*number_of_links_f64 - 1.0)*(sum_1/sum_0).powi(2) - (0.5*number_of_links_f64 - 1.5)*sum_2/sum_0) - nondimensional_potential_distance.powi(-2) - ((0.5*number_of_links_f64 - 1.0)*sum_1/sum_0 + nondimensional_potential_distance.powi(-1)).powi(2))
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
