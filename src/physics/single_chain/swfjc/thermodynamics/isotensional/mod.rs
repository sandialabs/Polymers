mod test;
pub mod legendre;
use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::ZERO;
pub struct SWFJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub well_width: f64,
    number_of_links_f64: f64,
    pub nondimensional_well_parameter: f64,
    pub legendre: self::legendre::SWFJC
}
impl SWFJC
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, well_width: f64) -> Self
    {
        SWFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            well_width,
            number_of_links_f64: number_of_links as f64,
            nondimensional_well_parameter: 1.0 + well_width/link_length,
            legendre: self::legendre::SWFJC::init(number_of_links, link_length, hinge_mass, well_width)
        }
    }
    pub fn end_to_end_length(&self, force: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = force*self.link_length/BOLTZMANN_CONSTANT/temperature;
        self.number_of_links_f64*self.link_length*((self.nondimensional_well_parameter.powi(2)*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.sinh())/(self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()) - 3.0/nondimensional_force)
    }
    pub fn end_to_end_length_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = force*self.link_length/BOLTZMANN_CONSTANT/temperature;
        self.link_length*((self.nondimensional_well_parameter.powi(2)*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.sinh())/(self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()) - 3.0/nondimensional_force)
    }
    pub fn nondimensional_end_to_end_length(&self, nondimensional_force: &f64) -> f64
    {
        self.number_of_links_f64*((self.nondimensional_well_parameter.powi(2)*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.sinh())/(self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()) - 3.0/nondimensional_force)
    }
    pub fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: &f64) -> f64
    {
        (self.nondimensional_well_parameter.powi(2)*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.sinh())/(self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()) - 3.0/nondimensional_force
    }
    pub fn gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = force*self.link_length/BOLTZMANN_CONSTANT/temperature;
        self.number_of_links_f64*BOLTZMANN_CONSTANT*temperature*(3.0*nondimensional_force.ln() - (self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()).ln()) - self.number_of_links_f64*BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = force*self.link_length/BOLTZMANN_CONSTANT/temperature;
        BOLTZMANN_CONSTANT*temperature*(3.0*nondimensional_force.ln() - (self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()).ln()) - BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn relative_gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        self.gibbs_free_energy(force, temperature) - self.gibbs_free_energy(&(ZERO*BOLTZMANN_CONSTANT*temperature/self.link_length), temperature)
    }
    pub fn relative_gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        self.gibbs_free_energy_per_link(force, temperature) - self.gibbs_free_energy_per_link(&(ZERO*BOLTZMANN_CONSTANT*temperature/self.link_length), temperature)
    }
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        self.number_of_links_f64*(3.0*nondimensional_force.ln() - (self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()).ln()) - self.number_of_links_f64*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        3.0*nondimensional_force.ln() - (self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()).ln() - (8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: &f64) -> f64
    {
        self.nondimensional_gibbs_free_energy(nondimensional_force, &300.0) - self.nondimensional_gibbs_free_energy(&ZERO, &300.0)
    }
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: &f64) -> f64
    {
        self.nondimensional_gibbs_free_energy_per_link(nondimensional_force, &300.0) - self.nondimensional_gibbs_free_energy_per_link(&ZERO, &300.0)
    }
}
