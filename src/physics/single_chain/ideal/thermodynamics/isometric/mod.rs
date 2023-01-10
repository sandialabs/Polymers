mod test;
use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
pub struct Ideal
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    number_of_links_f64: f64,
    contour_length: f64
}
impl Ideal
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        Ideal
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length
        }
    }
    pub fn force(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        3.0*end_to_end_length*BOLTZMANN_CONSTANT*temperature/self.number_of_links_f64/self.link_length.powi(2)
    }
    pub fn nondimensional_force(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        3.0*nondimensional_end_to_end_length_per_link
    }
    pub fn helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        self.helmholtz_free_energy_per_link(end_to_end_length, temperature)*self.number_of_links_f64
    }
    pub fn helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy_per_link(&(end_to_end_length/self.contour_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    pub fn relative_helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        self.relative_helmholtz_free_energy_per_link(end_to_end_length, temperature)*self.number_of_links_f64
    }
    pub fn relative_helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(&(end_to_end_length/self.contour_length))*BOLTZMANN_CONSTANT*temperature
    }
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        1.5*self.number_of_links_f64*nondimensional_end_to_end_length_per_link.powi(2) - self.number_of_links_f64*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        1.5*nondimensional_end_to_end_length_per_link.powi(2) - (8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        1.5*self.number_of_links_f64*nondimensional_end_to_end_length_per_link.powi(2)
    }
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        1.5*nondimensional_end_to_end_length_per_link.powi(2)
    }
    pub fn equilibrium_distribution(&self, end_to_end_length: &f64) -> f64
    {
        (1.5/PI/self.number_of_links_f64/self.link_length.powi(2)).powf(1.5)*(-1.5*(end_to_end_length/self.link_length).powi(2)/self.number_of_links_f64).exp()
    }
    pub fn nondimensional_equilibrium_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        (1.5/PI*self.number_of_links_f64).powf(1.5)*(-1.5*nondimensional_end_to_end_length_per_link.powi(2)*self.number_of_links_f64).exp()
    }
    pub fn equilibrium_radial_distribution(&self, end_to_end_length: &f64) -> f64
    {
        self.nondimensional_equilibrium_radial_distribution(&(end_to_end_length/self.contour_length))/self.contour_length
    }
    pub fn nondimensional_equilibrium_radial_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        4.0*PI*nondimensional_end_to_end_length_per_link.powi(2)*(1.5/PI*self.number_of_links_f64).powf(1.5)*(-1.5*nondimensional_end_to_end_length_per_link.powi(2)*self.number_of_links_f64).exp()
    }
}
