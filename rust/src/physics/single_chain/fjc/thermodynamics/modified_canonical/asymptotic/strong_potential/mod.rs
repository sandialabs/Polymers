pub mod test;
use std::f64::consts::PI;
use crate::math::
{
    ln,
    binomial,
    factorial
};
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::fjc::thermodynamics::ModifiedCanonicalAsymptoticStrongPotential;
use crate::physics::single_chain::fjc::ZERO;
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub number_of_links_f64: f64,
    pub contour_length: f64
}
impl ModifiedCanonicalAsymptoticStrongPotential for FJC
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
        BOLTZMANN_CONSTANT*temperature/self.link_length*self.nondimensional_force(&(potential_distance/self.contour_length), &(potential_stiffness*(self.contour_length).powf(2.0)/BOLTZMANN_CONSTANT/temperature))
    }
    fn nondimensional_force(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        1e-6*nondimensional_potential_distance
    }
    fn helmholtz_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        BOLTZMANN_CONSTANT*temperature*self.nondimensional_helmholtz_free_energy(&(potential_distance/self.contour_length), &(potential_stiffness*(self.contour_length).powf(2.0)/BOLTZMANN_CONSTANT/temperature), temperature)
    }
    fn helmholtz_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        BOLTZMANN_CONSTANT*temperature*self.nondimensional_helmholtz_free_energy_per_link(&(potential_distance/self.contour_length), &(potential_stiffness*(self.contour_length).powf(2.0)/BOLTZMANN_CONSTANT/temperature), temperature)
    }
    fn relative_helmholtz_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        BOLTZMANN_CONSTANT*temperature*self.nondimensional_relative_helmholtz_free_energy(&(potential_distance/self.contour_length), &(potential_stiffness*(self.contour_length).powf(2.0)/BOLTZMANN_CONSTANT/temperature))
    }
    fn relative_helmholtz_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        BOLTZMANN_CONSTANT*temperature*self.nondimensional_relative_helmholtz_free_energy_per_link(&(potential_distance/self.contour_length), &(potential_stiffness*(self.contour_length).powf(2.0)/BOLTZMANN_CONSTANT/temperature))
    }
    fn nondimensional_helmholtz_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let n = self.number_of_links as u128;
        let p = self.number_of_links_f64 - 2.0;
        let m = -*nondimensional_potential_distance*0.5 + 0.5;
        let k = (self.number_of_links_f64*m).ceil() as u128;
        let sum_0: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(binomial(&n, s) as f64)*(m - (*s as f64)/self.number_of_links_f64).powf(p)).sum();
        let sum_1: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(binomial(&n, s) as f64)*(m - (*s as f64)/self.number_of_links_f64).powf(p - 1.0)).sum();
        let sum_2: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(binomial(&n, s) as f64)*(m - (*s as f64)/self.number_of_links_f64).powf(p - 2.0)).sum();
        -ln(&(0.125/PI/nondimensional_potential_distance*(n.pow(n as u32) as f64)/(factorial(&(n - 2)) as f64)*sum_0/self.contour_length.powf(3.0))) - self.number_of_links_f64*ln(&(8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)))

        - 1.5*ln(&(2.0*PI/nondimensional_potential_stiffness)) - 3.0*ln(&(self.contour_length))
        
        + 0.5/nondimensional_potential_stiffness*(

            (0.5*self.number_of_links_f64 - 1.0)*(

                (0.5*self.number_of_links_f64 - 1.0)*(sum_1/sum_0).powf(2.0) - (0.5*self.number_of_links_f64 - 1.5)*sum_2/sum_0

            )
            
            - nondimensional_potential_distance.powf(-2.0)
            
            - (
                
                (0.5*self.number_of_links_f64 - 1.0)*sum_1/sum_0 + nondimensional_potential_distance.powf(-1.0)
            
            ).powf(2.0)
        )
    }
    fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness, temperature)/self.number_of_links_f64
    }
    fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness, &300.0) - self.nondimensional_helmholtz_free_energy(&ZERO, nondimensional_potential_stiffness, &300.0)
    }
    fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness)/self.number_of_links_f64
    }
}