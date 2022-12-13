pub mod test;
pub mod legendre;
use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::fjc::ZERO;
pub struct SWFJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub well_width: f64,
    pub number_of_links_f64: f64,
    pub contour_length: f64,
    pub nondimensional_well_parameter: f64,
    pub legendre: self::legendre::SWFJC
}
use super::Isotensional;
impl Isotensional for SWFJC
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, well_width: f64) -> SWFJC
    {
        SWFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            well_width,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length,
            nondimensional_well_parameter: 1.0 + well_width/link_length,
            legendre: self::legendre::SWFJC::init(number_of_links, link_length, hinge_mass, well_width)
        }
    }
    fn end_to_end_length(&self, force: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = force*self.link_length/BOLTZMANN_CONSTANT/temperature;
        self.number_of_links_f64*self.link_length*((self.nondimensional_well_parameter.powf(2.0)*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.sinh())/(self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()) - 3.0/nondimensional_force)
    }
    fn end_to_end_length_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = force*self.link_length/BOLTZMANN_CONSTANT/temperature;
        self.link_length*((self.nondimensional_well_parameter.powf(2.0)*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.sinh())/(self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()) - 3.0/nondimensional_force)
    }
    fn nondimensional_end_to_end_length(&self, nondimensional_force: &f64) -> f64
    {
        self.number_of_links_f64*((self.nondimensional_well_parameter.powf(2.0)*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.sinh())/(self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()) - 3.0/nondimensional_force)
    }
    fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: &f64) -> f64
    {
        (self.nondimensional_well_parameter.powf(2.0)*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.sinh())/(self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()) - 3.0/nondimensional_force
    }
    fn gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = force*self.link_length/BOLTZMANN_CONSTANT/temperature;
        self.number_of_links_f64*BOLTZMANN_CONSTANT*temperature*(3.0*nondimensional_force.ln() - (self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()).ln()) - self.number_of_links_f64*BOLTZMANN_CONSTANT*temperature*(8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)).ln()
    }
    fn gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = force*self.link_length/BOLTZMANN_CONSTANT/temperature;
        BOLTZMANN_CONSTANT*temperature*(3.0*nondimensional_force.ln() - (self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()).ln()) - BOLTZMANN_CONSTANT*temperature*(8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)).ln()
    }
    fn relative_gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        self.gibbs_free_energy(force, temperature) - self.gibbs_free_energy(&(ZERO*BOLTZMANN_CONSTANT*temperature/self.link_length), temperature)
    }
    fn relative_gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        self.gibbs_free_energy_per_link(force, temperature) - self.gibbs_free_energy_per_link(&(ZERO*BOLTZMANN_CONSTANT*temperature/self.link_length), temperature)
    }
    fn nondimensional_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        self.number_of_links_f64*(3.0*nondimensional_force.ln() - (self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()).ln()) - self.number_of_links_f64*(8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)).ln()
    }
    fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        3.0*nondimensional_force.ln() - (self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()).ln() - (8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)).ln()
    }
    fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: &f64) -> f64
    {
        self.nondimensional_gibbs_free_energy(nondimensional_force, &300.0) - self.nondimensional_gibbs_free_energy(&ZERO, &300.0)
    }
    fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: &f64) -> f64
    {
        self.nondimensional_gibbs_free_energy_per_link(nondimensional_force, &300.0) - self.nondimensional_gibbs_free_energy_per_link(&ZERO, &300.0)
    }
}
pub trait Legendre
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, well_width: f64) -> Self;
    fn helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64;
    fn helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64;
    fn relative_helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64;
    fn relative_helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64;
    fn nondimensional_helmholtz_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
    fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
    fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_force: &f64) -> f64;
    fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64) -> f64;
}
