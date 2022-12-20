pub mod test;
pub mod legendre;
use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::fjc::ZERO;
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub number_of_links_f64: f64,
    pub contour_length: f64,
    pub legendre: legendre::FJC
}
use super::Isometric;
impl Isometric for FJC
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> FJC
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length,
            legendre: legendre::FJC::init(number_of_links, link_length, hinge_mass)
        }
    }
    fn force(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_force(&(end_to_end_length/self.contour_length))*BOLTZMANN_CONSTANT*temperature/self.link_length
    }
    fn nondimensional_force(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        let n = self.number_of_links as u128;
        let p = self.number_of_links_f64 - 2.0;
        let m = -*nondimensional_end_to_end_length_per_link*0.5 + 0.5;
        let k = (self.number_of_links_f64*m).ceil() as u128;
        let sum_0: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powf(p)).sum();
        let sum_1: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powf(p - 1.0)).sum();
        (1.0/nondimensional_end_to_end_length_per_link + (0.5*self.number_of_links_f64 - 1.0)*sum_1/sum_0)/self.number_of_links_f64
    }
    fn helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        self.helmholtz_free_energy_per_link(end_to_end_length, temperature)*self.number_of_links_f64
    }
    fn helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy_per_link(&(end_to_end_length/self.contour_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    fn relative_helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        self.relative_helmholtz_free_energy_per_link(end_to_end_length, temperature)*self.number_of_links_f64
    }
    fn relative_helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(&(end_to_end_length/self.contour_length))*BOLTZMANN_CONSTANT*temperature
    }
    fn nondimensional_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        -(self.equilibrium_distribution(&(nondimensional_end_to_end_length_per_link*self.contour_length))).ln() - (self.number_of_links_f64 - 1.0)*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(nondimensional_end_to_end_length_per_link, temperature)/self.number_of_links_f64
    }
    fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        (self.nondimensional_equilibrium_distribution(&ZERO)/self.nondimensional_equilibrium_distribution(nondimensional_end_to_end_length_per_link)).ln()
    }
    fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy(nondimensional_end_to_end_length_per_link)/self.number_of_links_f64
    }
    fn equilibrium_distribution(&self, end_to_end_length: &f64) -> f64
    {
        self.nondimensional_equilibrium_distribution(&(end_to_end_length/self.contour_length))/self.contour_length.powi(3)
    }
    fn nondimensional_equilibrium_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        let n = self.number_of_links as u128;
        let p = self.number_of_links_f64 - 2.0;
        let m = -*nondimensional_end_to_end_length_per_link*0.5 + 0.5;
        let k = (self.number_of_links_f64*m).ceil() as u128;
        let sum: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powf(p)).sum();
        0.125/PI/nondimensional_end_to_end_length_per_link*(n.pow(n as u32) as f64)/((1..=n-2).product::<u128>() as f64)*sum
    }
    fn equilibrium_radial_distribution(&self, end_to_end_length: &f64) -> f64
    {
        self.nondimensional_equilibrium_radial_distribution(&(end_to_end_length/self.contour_length))/self.contour_length
    }
    fn nondimensional_equilibrium_radial_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        let n = self.number_of_links as u128;
        let p = self.number_of_links_f64 - 2.0;
        let m = -*nondimensional_end_to_end_length_per_link*0.5 + 0.5;
        let k = (self.number_of_links_f64*m).ceil() as u128;
        let sum: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powf(p)).sum();
        0.5*nondimensional_end_to_end_length_per_link*(n.pow(n as u32) as f64)/((1..=n-2).product::<u128>() as f64)*sum
    }
}
pub trait Legendre
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self;
    fn force(&self, end_to_end_length: &f64, temperature: &f64) -> f64;
    fn nondimensional_force(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64;
    fn helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64;
    fn helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64;
    fn relative_helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64;
    fn relative_helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64;
    fn nondimensional_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64;
    fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64;
    fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64;
    fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64;
    fn equilibrium_distribution(&self, end_to_end_length: &f64) -> f64;
    fn nondimensional_equilibrium_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64;
    fn equilibrium_radial_distribution(&self, end_to_end_length: &f64) -> f64;
    fn nondimensional_equilibrium_radial_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64;
    fn gibbs_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64;
    fn gibbs_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64;
    fn relative_gibbs_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64;
    fn relative_gibbs_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64;
    fn nondimensional_gibbs_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64;
    fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64;
    fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64;
    fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64;
}