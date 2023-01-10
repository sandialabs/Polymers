mod test;
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
    pub nondimensional_well_parameter: f64
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

        }
    }
    pub fn helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(&(*force/BOLTZMANN_CONSTANT/temperature*self.link_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    pub fn helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy_per_link(&(*force/BOLTZMANN_CONSTANT/temperature*self.link_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    pub fn relative_helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy(&(*force/BOLTZMANN_CONSTANT/temperature*self.link_length))*BOLTZMANN_CONSTANT*temperature
    }
    pub fn relative_helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(&(*force/BOLTZMANN_CONSTANT/temperature*self.link_length))*BOLTZMANN_CONSTANT*temperature
    }
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        self.number_of_links_f64*self.nondimensional_helmholtz_free_energy_per_link(nondimensional_force, temperature)
    }
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
    (self.nondimensional_well_parameter.powi(2)*nondimensional_force.powi(2)*(self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force.powi(2)*nondimensional_force.sinh())/(self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()) - 3.0 + 3.0*nondimensional_force.ln() - (self.nondimensional_well_parameter*nondimensional_force*(self.nondimensional_well_parameter*nondimensional_force).cosh() - (self.nondimensional_well_parameter*nondimensional_force).sinh() - nondimensional_force*nondimensional_force.cosh() + nondimensional_force.sinh()).ln() - (8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_force: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(nondimensional_force, &300.0) - self.nondimensional_helmholtz_free_energy(&ZERO, &300.0)
    }
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy_per_link(nondimensional_force, &300.0) - self.nondimensional_helmholtz_free_energy_per_link(&ZERO, &300.0)
    }
}
