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

/// The structure of the Morse-FJC model thermodynamics in the isotensional ensemble approximated using a Legendre transformation.
pub struct MORSEFJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The stiffness of each link in the chain in units of J/(mol⋅nm^2).
    pub link_stiffness: f64,

    /// The energy of each link in the chain in units of J/mol.
    pub link_energy: f64
}

/// The Helmholtz free energy as a function of the applied force and temperature, parameterized by the number of links, link length, hinge mass, link stiffness, and link stiffness.
pub fn helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, link_stiffness: &f64, link_energy: &f64, force: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(force/BOLTZMANN_CONSTANT/temperature*link_length), temperature)
}

/// The Helmholtz free energy per link as a function of the applied force and temperature, parameterized by the link length, hinge mass, link stiffness, and link stiffness.
pub fn helmholtz_free_energy_per_link(link_length: &f64, hinge_mass: &f64, link_stiffness: &f64, link_energy: &f64, force: &f64, temperature: &f64) -> f64
{
    BOLTZMANN_CONSTANT*temperature*nondimensional_helmholtz_free_energy_per_link(link_length, hinge_mass, &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(force/BOLTZMANN_CONSTANT/temperature*link_length), temperature)
}

/// The relative Helmholtz free energy as a function of the applied force and temperature, parameterized by the number of links, link length, link stiffness, and link stiffness.
pub fn relative_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, link_stiffness: &f64, link_energy: &f64, force: &f64, temperature: &f64) -> f64
{
    helmholtz_free_energy(number_of_links, link_length, &1.0, link_stiffness, link_energy, force, temperature) - helmholtz_free_energy(number_of_links, link_length, &1.0, link_stiffness, link_energy, &(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), temperature)
}

/// The relative Helmholtz free energy per link as a function of the applied force and temperature, parameterized by the link length link stiffness, and link stiffness.
pub fn relative_helmholtz_free_energy_per_link(link_length: &f64, link_stiffness: &f64, link_energy: &f64, force: &f64, temperature: &f64) -> f64
{
    helmholtz_free_energy_per_link(link_length, &1.0, link_stiffness, link_energy, force, temperature) - helmholtz_free_energy_per_link(link_length, &1.0, link_stiffness, link_energy, &(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), temperature)
}

/// The nondimensional Helmholtz free energy as a function of the applied nondimensional force and temperature, parameterized by the number of links, link length, hinge mass, nondimensional link stiffness, and nondimensional link energy.
pub fn nondimensional_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_link_stiffness: &f64, nondimensional_link_energy: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    (*number_of_links as f64)*nondimensional_helmholtz_free_energy_per_link(link_length, hinge_mass, nondimensional_link_stiffness, nondimensional_link_energy, nondimensional_force, temperature)
}

/// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional force and temperature, parameterized by the link length, hinge mass, nondimensional link stiffness, and nondimensional link energy.
pub fn nondimensional_helmholtz_free_energy_per_link(link_length: &f64, hinge_mass: &f64, nondimensional_link_stiffness: &f64, nondimensional_link_energy: &f64, nondimensional_force: &f64, temperature: &f64) -> f64
{
    let nondimensional_morse_parameter = (nondimensional_link_stiffness/nondimensional_link_energy/2.0).sqrt();
    let nondimensional_link_stretch_max = 1.0 + 2.0_f64.ln()/nondimensional_morse_parameter;
    let rescaled_partition_function_integrand = |nondimensional_link_stretch: f64|
    {
        let exponent_1 = nondimensional_force*nondimensional_link_stretch - nondimensional_link_energy*(1.0 - (-nondimensional_morse_parameter*(nondimensional_link_stretch - 1.0)).exp()).powi(2) + nondimensional_link_stretch.ln() - nondimensional_force.ln();
        let exponent_2 = exponent_1 - 2.0*nondimensional_force*nondimensional_link_stretch;
        exponent_1.exp() - exponent_2.exp()
    };
    let dx = (ONE*nondimensional_link_stretch_max - ZERO)/(POINTS as f64);
    let rescaled_partition_function = (0..=POINTS-1).collect::<Vec::<u128>>().iter().map(|index| rescaled_partition_function_integrand(ZERO + (0.5 + *index as f64)*dx)).sum::<f64>();
    let nondimensional_end_to_end_length_per_link_integrand = |nondimensional_link_stretch: f64|
    {
        let exponent_1 = nondimensional_force*nondimensional_link_stretch - nondimensional_link_energy*(1.0 - (-nondimensional_morse_parameter*(nondimensional_link_stretch - 1.0)).exp()).powi(2) + 2.0*nondimensional_link_stretch.ln() - nondimensional_force.ln();
        let exponent_2 = exponent_1 - 2.0*nondimensional_force*nondimensional_link_stretch;
        let exponent_3 = exponent_1 - nondimensional_link_stretch.ln() - nondimensional_force.ln();
        let exponent_4 = exponent_3 - 2.0*nondimensional_force*nondimensional_link_stretch;
        (exponent_1.exp() + exponent_2.exp() - exponent_3.exp() + exponent_4.exp())/rescaled_partition_function
    };
    let nondimensional_end_to_end_length_per_link = (0..=POINTS-1).collect::<Vec::<u128>>().iter().map(|index| nondimensional_end_to_end_length_per_link_integrand(ZERO + (0.5 + *index as f64)*dx)).sum::<f64>();
    -rescaled_partition_function.ln() - (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln() + nondimensional_force*nondimensional_end_to_end_length_per_link
}

/// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional force, parameterized by the number of links and nondimensional link stiffness and nondimensional link energy.
pub fn nondimensional_relative_helmholtz_free_energy(number_of_links: &u8, nondimensional_link_stiffness: &f64, nondimensional_link_energy: &f64, nondimensional_force: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, &1.0, &1.0, nondimensional_link_stiffness, nondimensional_link_energy, nondimensional_force, &300.0) - nondimensional_helmholtz_free_energy(number_of_links, &1.0, &1.0, nondimensional_link_stiffness, nondimensional_link_energy, &ZERO, &300.0)
}

/// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional force, parameterized by the nondimensional link stiffness and nondimensional link energy.
pub fn nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_link_stiffness: &f64, nondimensional_link_energy: &f64, nondimensional_force: &f64) -> f64
{
    nondimensional_helmholtz_free_energy_per_link(&1.0, &1.0, nondimensional_link_stiffness, nondimensional_link_energy, nondimensional_force, &300.0) - nondimensional_helmholtz_free_energy_per_link(&1.0, &1.0, nondimensional_link_stiffness, nondimensional_link_energy, &ZERO, &300.0)
}

/// The implemented functionality of the Morse-FJC model thermodynamics in the isotensional ensemble approximated using a Legendre transformation.
impl MORSEFJC
{
    /// Initializes and returns an instance of the Morse-FJC model thermodynamics in the isotensional ensemble approximated using a Legendre transformation.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64, link_energy: f64) -> Self
    {
        MORSEFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            link_energy
        }
    }
    /// The Helmholtz free energy as a function of the applied force and temperature.
    pub fn helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.link_stiffness, &self.link_energy, force, temperature)
    }
    /// The Helmholtz free energy per link as a function of the applied force and temperature.
    pub fn helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy_per_link(&self.link_length, &self.hinge_mass, &self.link_stiffness, &self.link_energy, force, temperature)
    }
    /// The relative Helmholtz free energy as a function of the applied force and temperature.
    pub fn relative_helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.link_stiffness, &self.link_energy, force, temperature)
    }
    /// The relative Helmholtz free energy per link as a function of the applied force and temperature.
    pub fn relative_helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy_per_link(&self.link_length, &self.link_stiffness, &self.link_energy, force, temperature)
    }
    /// The nondimensional Helmholtz free energy as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(self.link_energy/BOLTZMANN_CONSTANT/temperature), nondimensional_force, temperature)
    }
    /// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional force and temperature.
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy_per_link(&self.link_length, &self.hinge_mass, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(self.link_energy/BOLTZMANN_CONSTANT/temperature), nondimensional_force, temperature)
    }
    /// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional force.
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy(&self.number_of_links, &(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(self.link_energy/BOLTZMANN_CONSTANT/temperature), nondimensional_force)
    }
    /// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional force.
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy_per_link(&(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(self.link_energy/BOLTZMANN_CONSTANT/temperature), nondimensional_force)
    }
}
