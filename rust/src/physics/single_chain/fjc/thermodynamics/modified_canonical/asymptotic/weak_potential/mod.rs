pub mod test;
use std::f64::consts::PI;
use crate::math::
{
    ln,
    ln_sinhc,
    langevin
};
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::fjc::thermodynamics::ModifiedCanonicalAsymptoticWeakPotential;
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub number_of_links_f64: f64,
    pub contour_length: f64
}
impl ModifiedCanonicalAsymptoticWeakPotential for FJC
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> FJC
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length
        }
    }
    fn force(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        potential_stiffness*potential_distance
    }
    fn nondimensional_force(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64
    }
    fn end_to_end_length(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = potential_stiffness*potential_distance*self.link_length/BOLTZMANN_CONSTANT/temperature;
        self.contour_length*langevin(&nondimensional_force) - potential_stiffness*self.contour_length.powf(2.0)*self.link_length/BOLTZMANN_CONSTANT/temperature*(langevin(&nondimensional_force)*(nondimensional_force.powf(-2.0) - (nondimensional_force.sinh()).powf(-2.0)) + ((nondimensional_force.sinh()).powf(-2.0)/nondimensional_force.tanh() - nondimensional_force.powf(-3.0))/self.number_of_links_f64)
    }
    fn end_to_end_length_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = potential_stiffness*potential_distance*self.link_length/BOLTZMANN_CONSTANT/temperature;
        self.link_length*langevin(&nondimensional_force) - potential_stiffness*self.number_of_links_f64*self.link_length.powf(3.0)/BOLTZMANN_CONSTANT/temperature*(langevin(&nondimensional_force)*(nondimensional_force.powf(-2.0) - (nondimensional_force.sinh()).powf(-2.0)) + ((nondimensional_force.sinh()).powf(-2.0)/nondimensional_force.tanh() - nondimensional_force.powf(-3.0))/self.number_of_links_f64)
    }
    fn nondimensional_end_to_end_length(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        let nondimensional_force = nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64;
        self.number_of_links_f64*langevin(&nondimensional_force) - nondimensional_potential_stiffness*(langevin(&nondimensional_force)*(nondimensional_force.powf(-2.0) - (nondimensional_force.sinh()).powf(-2.0)) + ((nondimensional_force.sinh()).powf(-2.0)/nondimensional_force.tanh() - nondimensional_force.powf(-3.0))/self.number_of_links_f64)
    }
    fn nondimensional_end_to_end_length_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        let nondimensional_force = nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64;
        langevin(&nondimensional_force) - nondimensional_potential_stiffness*(langevin(&nondimensional_force)*(nondimensional_force.powf(-2.0) - (nondimensional_force.sinh()).powf(-2.0)) + ((nondimensional_force.sinh()).powf(-2.0)/nondimensional_force.tanh() - nondimensional_force.powf(-3.0))/self.number_of_links_f64)/self.number_of_links_f64
    }
    fn gibbs_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.relative_gibbs_free_energy(potential_distance, potential_stiffness, temperature) - self.number_of_links_f64*BOLTZMANN_CONSTANT*temperature*ln(&(8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)))
    }
    fn gibbs_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.relative_gibbs_free_energy_per_link(potential_distance, potential_stiffness, temperature) - BOLTZMANN_CONSTANT*temperature*ln(&(8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)))
    }
    fn relative_gibbs_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = potential_stiffness*potential_distance*self.link_length/BOLTZMANN_CONSTANT/temperature;
        -self.number_of_links_f64*BOLTZMANN_CONSTANT*temperature*ln_sinhc(&nondimensional_force) + 0.5*potential_stiffness*self.contour_length.powf(2.0)*(langevin(&nondimensional_force).powf(2.0) + (nondimensional_force.powf(-2.0) - (nondimensional_force.sinh()).powf(-2.0) - 1.0/3.0)/self.number_of_links_f64)
    }
    fn relative_gibbs_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let nondimensional_force = potential_stiffness*potential_distance*self.link_length/BOLTZMANN_CONSTANT/temperature;
        -BOLTZMANN_CONSTANT*temperature*ln_sinhc(&nondimensional_force) + 0.5*potential_stiffness*self.contour_length.powf(2.0)*(langevin(&nondimensional_force).powf(2.0) + (nondimensional_force.powf(-2.0) - (nondimensional_force.sinh()).powf(-2.0) - 1.0/3.0)/self.number_of_links_f64)/self.number_of_links_f64
    }
    fn nondimensional_gibbs_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_gibbs_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness) - self.number_of_links_f64*ln(&(8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)))
    }
    fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_gibbs_free_energy_per_link(nondimensional_potential_distance, nondimensional_potential_stiffness) - ln(&(8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)))
    }
    fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        let nondimensional_force = nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64;
        -self.number_of_links_f64*ln_sinhc(&nondimensional_force) + 0.5*nondimensional_potential_stiffness*(langevin(&nondimensional_force).powf(2.0) + (nondimensional_force.powf(-2.0) - (nondimensional_force.sinh()).powf(-2.0) - 1.0/3.0)/self.number_of_links_f64)
    }
    fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        let nondimensional_force = nondimensional_potential_stiffness*nondimensional_potential_distance/self.number_of_links_f64;
        -ln_sinhc(&nondimensional_force) + 0.5*nondimensional_potential_stiffness*(langevin(&nondimensional_force).powf(2.0) + (nondimensional_force.powf(-2.0) - (nondimensional_force.sinh()).powf(-2.0) - 1.0/3.0)/self.number_of_links_f64)/self.number_of_links_f64
    }
}