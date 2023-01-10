mod test;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT,
    single_chain::ZERO
};
use std::f64::consts::PI;
pub struct EFJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub link_stiffness: f64,
    number_of_links_f64: f64
}
impl EFJC
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        EFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            number_of_links_f64: number_of_links as f64
        }
    }
    pub fn helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(&(force/BOLTZMANN_CONSTANT/temperature*self.link_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    pub fn helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy_per_link(&(force/BOLTZMANN_CONSTANT/temperature*self.link_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    pub fn relative_helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        self.helmholtz_free_energy(force, temperature) - self.helmholtz_free_energy(&(ZERO*BOLTZMANN_CONSTANT*temperature/self.link_length), temperature)
    }
    pub fn relative_helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        self.helmholtz_free_energy_per_link(force, temperature) - self.helmholtz_free_energy_per_link(&(ZERO*BOLTZMANN_CONSTANT*temperature/self.link_length), temperature)
    }
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        self.number_of_links_f64*self.nondimensional_helmholtz_free_energy_per_link(nondimensional_force, temperature)
    }
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        -(nondimensional_force.sinh()/nondimensional_force).ln() - 0.5*nondimensional_force.powi(2)/(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature) - 0.5*(2.0*PI*BOLTZMANN_CONSTANT*temperature/self.link_stiffness).ln() - (8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln() + nondimensional_force/nondimensional_force.tanh() - 1.0 + nondimensional_force.powi(2)/(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature)
    }
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(nondimensional_force, temperature) - self.nondimensional_helmholtz_free_energy(&ZERO, temperature)
    }
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy_per_link(nondimensional_force, temperature) - self.nondimensional_helmholtz_free_energy_per_link(&ZERO, temperature)
    }
}
