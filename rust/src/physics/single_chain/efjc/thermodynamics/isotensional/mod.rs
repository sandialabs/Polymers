pub mod test;
pub mod asymptotic;
pub mod legendre;
use self::asymptotic::{
    Alternative as AsymptoticAlternative,
    Reduced as AsymptoticReduced,
    alternative::Legendre as AsymptoticAlternativeLegendre,
    reduced::Legendre as AsymptoticReducedLegendre
};
use std::f64::consts::PI;
use crate::math::erf;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::fjc::ZERO;
pub struct EFJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub link_stiffness: f64,
    pub number_of_links_f64: f64,
    pub contour_length: f64,
    pub asymptotic: self::asymptotic::EFJC,
    pub legendre: self::legendre::EFJC
}
use super::Isotensional;
impl Isotensional for EFJC
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> EFJC
    {
        EFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length,
            asymptotic: self::asymptotic::EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness),
            legendre: self::legendre::EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness)
        }
    }
    fn end_to_end_length(&self, force: &f64, temperature: &f64) -> f64
    {
        self.link_length*self.nondimensional_end_to_end_length(&(force*self.link_length/BOLTZMANN_CONSTANT/temperature), temperature)
    }
    fn end_to_end_length_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        self.link_length*self.nondimensional_end_to_end_length_per_link(&(force*self.link_length/BOLTZMANN_CONSTANT/temperature), temperature)
    }
    fn nondimensional_end_to_end_length(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        self.number_of_links_f64*self.nondimensional_end_to_end_length_per_link(nondimensional_force, temperature)
    }
    fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        let nondimensional_link_stiffness = self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature;
        let denominator = 4.0*nondimensional_force.sinh()*(1.0 + nondimensional_force/nondimensional_force.tanh()/nondimensional_link_stiffness);
        let fraction = ((nondimensional_force/nondimensional_link_stiffness + 1.0)*nondimensional_force.exp()*erf(&((nondimensional_force + nondimensional_link_stiffness)/(2.0*nondimensional_link_stiffness).sqrt())) - (nondimensional_force/nondimensional_link_stiffness - 1.0)/nondimensional_force.exp()*erf(&((nondimensional_force - nondimensional_link_stiffness)/(2.0*nondimensional_link_stiffness).sqrt())))/denominator;
        1.0/nondimensional_force.tanh() - 1.0/nondimensional_force + nondimensional_force/(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature)*(1.0 + (nondimensional_force.tanh() - 1.0/nondimensional_force.tanh() + 1.0/nondimensional_force)/(nondimensional_force.tanh() + nondimensional_force/(self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature))) + (nondimensional_force.exp()*((2.0/PI/nondimensional_link_stiffness).sqrt()*(nondimensional_force/nondimensional_link_stiffness + 1.0)*(-(nondimensional_force + nondimensional_link_stiffness).powi(2)/2.0/nondimensional_link_stiffness).exp() + (1.0 + (1.0 + nondimensional_force)/nondimensional_link_stiffness)) - 1.0/nondimensional_force.exp()*((2.0/PI/nondimensional_link_stiffness).sqrt()*(nondimensional_force/nondimensional_link_stiffness - 1.0)*(-(nondimensional_force - nondimensional_link_stiffness).powi(2)/2.0/nondimensional_link_stiffness).exp() + (1.0 + (1.0 - nondimensional_force)/nondimensional_link_stiffness)*erf(&((nondimensional_force - nondimensional_link_stiffness)/(2.0*nondimensional_link_stiffness).sqrt()))) - fraction*(4.0*(nondimensional_force.cosh()*(1.0 + (1.0 + nondimensional_force/nondimensional_force.tanh())/nondimensional_link_stiffness) - nondimensional_force/nondimensional_link_stiffness/nondimensional_force.sinh())))/denominator/(1.0 + fraction)
    }
    fn gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = force*self.link_length/BOLTZMANN_CONSTANT/temperature;
        let nondimensional_link_stiffness = self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature;
        self.number_of_links_f64*BOLTZMANN_CONSTANT*temperature*(-(nondimensional_force.sinh()/nondimensional_force).ln() - 0.5*nondimensional_force.powi(2)/nondimensional_link_stiffness - (1.0 + nondimensional_force/nondimensional_force.tanh()/nondimensional_link_stiffness).ln() - (0.5 + ((nondimensional_force/nondimensional_link_stiffness + 1.0)*nondimensional_force.exp()*erf(&((nondimensional_force + nondimensional_link_stiffness)/(2.0*nondimensional_link_stiffness).sqrt())) - (nondimensional_force/nondimensional_link_stiffness - 1.0)/nondimensional_force.exp()*erf(&((nondimensional_force - nondimensional_link_stiffness)/(2.0*nondimensional_link_stiffness).sqrt())))/(4.0*nondimensional_force.sinh()*(1.0 + nondimensional_force/nondimensional_force.tanh()/nondimensional_link_stiffness))).ln() - 0.5*(2.0*PI*BOLTZMANN_CONSTANT*temperature/self.link_stiffness).ln() - (8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())
    }
    fn gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = force*self.link_length/BOLTZMANN_CONSTANT/temperature;
        let nondimensional_link_stiffness = self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature;
        BOLTZMANN_CONSTANT*temperature*(-(nondimensional_force.sinh()/nondimensional_force).ln() - 0.5*nondimensional_force.powi(2)/nondimensional_link_stiffness - (1.0 + nondimensional_force/nondimensional_force.tanh()/nondimensional_link_stiffness).ln() - (0.5 + ((nondimensional_force/nondimensional_link_stiffness + 1.0)*nondimensional_force.exp()*erf(&((nondimensional_force + nondimensional_link_stiffness)/(2.0*nondimensional_link_stiffness).sqrt())) - (nondimensional_force/nondimensional_link_stiffness - 1.0)/nondimensional_force.exp()*erf(&((nondimensional_force - nondimensional_link_stiffness)/(2.0*nondimensional_link_stiffness).sqrt())))/(4.0*nondimensional_force.sinh()*(1.0 + nondimensional_force/nondimensional_force.tanh()/nondimensional_link_stiffness))).ln() - 0.5*(2.0*PI*BOLTZMANN_CONSTANT*temperature/self.link_stiffness).ln() - (8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())
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
        let nondimensional_link_stiffness = self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature;
        self.number_of_links_f64*(-(nondimensional_force.sinh()/nondimensional_force).ln() - 0.5*nondimensional_force.powi(2)/nondimensional_link_stiffness - (1.0 + nondimensional_force/nondimensional_force.tanh()/nondimensional_link_stiffness).ln() - (0.5 + ((nondimensional_force/nondimensional_link_stiffness + 1.0)*nondimensional_force.exp()*erf(&((nondimensional_force + nondimensional_link_stiffness)/(2.0*nondimensional_link_stiffness).sqrt())) - (nondimensional_force/nondimensional_link_stiffness - 1.0)/nondimensional_force.exp()*erf(&((nondimensional_force - nondimensional_link_stiffness)/(2.0*nondimensional_link_stiffness).sqrt())))/(4.0*nondimensional_force.sinh()*(1.0 + nondimensional_force/nondimensional_force.tanh()/nondimensional_link_stiffness))).ln() - 0.5*(2.0*PI*BOLTZMANN_CONSTANT*temperature/self.link_stiffness).ln() - (8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())
    }
    fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        let nondimensional_link_stiffness = self.link_stiffness*self.link_length.powi(2)/BOLTZMANN_CONSTANT/temperature;
        -(nondimensional_force.sinh()/nondimensional_force).ln() - 0.5*nondimensional_force.powi(2)/nondimensional_link_stiffness - (1.0 + nondimensional_force/nondimensional_force.tanh()/nondimensional_link_stiffness).ln() - (0.5 + ((nondimensional_force/nondimensional_link_stiffness + 1.0)*nondimensional_force.exp()*erf(&((nondimensional_force + nondimensional_link_stiffness)/(2.0*nondimensional_link_stiffness).sqrt())) - (nondimensional_force/nondimensional_link_stiffness - 1.0)/nondimensional_force.exp()*erf(&((nondimensional_force - nondimensional_link_stiffness)/(2.0*nondimensional_link_stiffness).sqrt())))/(4.0*nondimensional_force.sinh()*(1.0 + nondimensional_force/nondimensional_force.tanh()/nondimensional_link_stiffness))).ln() - 0.5*(2.0*PI*BOLTZMANN_CONSTANT*temperature/self.link_stiffness).ln() - (8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_gibbs_free_energy(nondimensional_force, temperature) - self.nondimensional_gibbs_free_energy(&ZERO, temperature)
    }
    fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_gibbs_free_energy_per_link(nondimensional_force, temperature) - self.nondimensional_gibbs_free_energy_per_link(&ZERO, temperature)
    }
}
pub trait Legendre
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self;
   fn helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64;
   fn helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64;
   fn relative_helmholtz_free_energy(&self, force: &f64, temperature: &f64) -> f64;
   fn relative_helmholtz_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64;
   fn nondimensional_helmholtz_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
   fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
   fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
   fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
}
pub trait Asymptotic
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self;
    fn end_to_end_length(&self, force: &f64, temperature: &f64) -> f64;
    fn end_to_end_length_per_link(&self, force: &f64, temperature: &f64) -> f64;
    fn nondimensional_end_to_end_length(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
    fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
    fn gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64;
    fn gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64;
    fn relative_gibbs_free_energy(&self, force: &f64, temperature: &f64) -> f64;
    fn relative_gibbs_free_energy_per_link(&self, force: &f64, temperature: &f64) -> f64;
    fn nondimensional_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
    fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
    fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
    fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: &f64) -> f64;
}
