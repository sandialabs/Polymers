#[cfg(feature = "extern")]
pub mod ex;

#[cfg(feature = "python")]
pub mod py;

mod test;

use super::
{
    nondimensional_force,
    nondimensional_helmholtz_free_energy
};
use crate::physics::BOLTZMANN_CONSTANT;
use crate::physics::single_chain::ZERO;

/// The structure of the thermodynamics of the WLC model in the isometric ensemble approximated using a Legendre transformation.
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

    nondimensional_persistance_length: f64
}

/// The Gibbs free energy as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, hinge mass, and persistance length.
pub fn gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, persistance_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_gibbs_free_energy(number_of_links, link_length, hinge_mass, &(persistance_length/contour_length), &(end_to_end_length/contour_length), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The Gibbs free energy per link as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, hinge mass, and persistance length.
pub fn gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, persistance_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_gibbs_free_energy_per_link(number_of_links, link_length, hinge_mass, &(persistance_length/contour_length), &(end_to_end_length/contour_length), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The relative Gibbs free energy as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, and persistance length.
pub fn relative_gibbs_free_energy(number_of_links: &u8, link_length: &f64, persistance_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_relative_gibbs_free_energy(number_of_links, &(persistance_length/contour_length), &(end_to_end_length/contour_length))*BOLTZMANN_CONSTANT*temperature
}

/// The relative Gibbs free energy per link as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, and persistance length.
pub fn relative_gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, persistance_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_relative_gibbs_free_energy_per_link(number_of_links, &(persistance_length/contour_length), &(end_to_end_length/contour_length))*BOLTZMANN_CONSTANT*temperature
}

/// The nondimensional Gibbs free energy as a function of the nondimensional end-to-end length per link and temperature, parameterized by the number of links, link length, hinge mass, and nondimensional persistance length.
pub fn nondimensional_gibbs_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_persistance_length: &f64, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, nondimensional_persistance_length, nondimensional_end_to_end_length_per_link, temperature) - nondimensional_force(number_of_links, nondimensional_persistance_length, nondimensional_end_to_end_length_per_link)*nondimensional_end_to_end_length_per_link*(*number_of_links as f64)
}

/// The nondimensional Gibbs free energy per link as a function of the applied nondimensional end-to-end length per link and temperature, parameterized by the number of links, link length, and nondimensional persistance length.
pub fn nondimensional_gibbs_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_persistance_length: &f64, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
{
    nondimensional_gibbs_free_energy(number_of_links, link_length, hinge_mass, nondimensional_persistance_length, nondimensional_end_to_end_length_per_link, temperature)/(*number_of_links as f64)
}

/// The nondimensional relative Gibbs free energy as a function of the nondimensional end-to-end length per link, parameterized by the number of links and nondimensional persistance length.
pub fn nondimensional_relative_gibbs_free_energy(number_of_links: &u8, nondimensional_persistance_length: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    nondimensional_gibbs_free_energy(number_of_links, &1.0, &1.0, nondimensional_persistance_length, nondimensional_end_to_end_length_per_link, &300.0) - nondimensional_gibbs_free_energy(number_of_links, &1.0, &1.0, nondimensional_persistance_length, &ZERO, &300.0)
}

/// The nondimensional relative Gibbs free energy as a function of the nondimensional end-to-end length per link, parameterized by the number of links and nondimensional persistance length.
pub fn nondimensional_relative_gibbs_free_energy_per_link(number_of_links: &u8, nondimensional_persistance_length: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    nondimensional_gibbs_free_energy_per_link(number_of_links, &1.0, &1.0, nondimensional_persistance_length, nondimensional_end_to_end_length_per_link, &300.0) - nondimensional_gibbs_free_energy_per_link(number_of_links, &1.0, &1.0, nondimensional_persistance_length, &ZERO, &300.0)
}

/// The implemented functionality of the thermodynamics of the WLC model in the isometric ensemble approximated using a Legendre transformation.
impl WLC
{
    /// Initializes and returns an instance of the thermodynamics of the WLC model in the isometric ensemble approximated using a Legendre transformation.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, persistance_length: f64) -> Self
    {
        WLC
        {
            hinge_mass,
            link_length,
            number_of_links,
            persistance_length,
            nondimensional_persistance_length: persistance_length/(number_of_links as f64)/link_length
        }
    }
    /// The Gibbs free energy as a function of the applied end-to-end length and temperature.
    pub fn gibbs_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.persistance_length, end_to_end_length, temperature)
    }
    /// The Gibbs free energy per link as a function of the applied end-to-end length and temperature.
    pub fn gibbs_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.persistance_length, end_to_end_length, temperature)
    }
    /// The relative Gibbs free energy as a function of the applied end-to-end length and temperature.
    pub fn relative_gibbs_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        relative_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.persistance_length, end_to_end_length, temperature)
    }
    /// The relative Gibbs free energy per link as a function of the applied end-to-end length and temperature.
    pub fn relative_gibbs_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        relative_gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.persistance_length, end_to_end_length, temperature)
    }
    /// The nondimensional Gibbs free energy as a function of the applied nondimensional end-to-end length per link and temperature.
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        nondimensional_gibbs_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.nondimensional_persistance_length, nondimensional_end_to_end_length_per_link, temperature)
    }
    /// The nondimensional Gibbs free energy per link as a function of the applied nondimensional end-to-end length per link and temperature.
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        nondimensional_gibbs_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.nondimensional_persistance_length, nondimensional_end_to_end_length_per_link, temperature)
    }
    /// The nondimensional relative Gibbs free energy as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_relative_gibbs_free_energy(&self.number_of_links, &self.nondimensional_persistance_length, nondimensional_end_to_end_length_per_link)
    }
    /// The nondimensional relative Gibbs free energy per link as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_relative_gibbs_free_energy_per_link(&self.number_of_links, &self.nondimensional_persistance_length, nondimensional_end_to_end_length_per_link)
    }
}
