pub mod test;
pub mod asymptotic;
use std::f64::consts::PI;
use crate::math::
{
    ln,
    factorial,
    integrate
};
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::fjc::thermodynamics::ModifiedCanonical;
pub static ONE: f64 = 1.0;
pub static ZERO: f64 = 1e-6;
pub static POINTS: u128 = 100;
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub number_of_links_f64: f64,
    pub contour_length: f64
}
impl ModifiedCanonical for FJC
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
    fn end_to_end_length(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        potential_distance - self.force(potential_distance, potential_stiffness, temperature)/potential_stiffness
    }
    fn end_to_end_length_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        (potential_distance - self.force(potential_distance, potential_stiffness, temperature)/potential_stiffness)/self.number_of_links_f64
    }
    fn nondimensional_end_to_end_length(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        self.number_of_links_f64*nondimensional_potential_distance - self.number_of_links_f64.powf(2.0)*self.nondimensional_force(nondimensional_potential_distance, nondimensional_potential_stiffness)/nondimensional_potential_stiffness
    }
    fn nondimensional_end_to_end_length_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        nondimensional_potential_distance - self.number_of_links_f64*self.nondimensional_force(nondimensional_potential_distance, nondimensional_potential_stiffness)/nondimensional_potential_stiffness
    }
    fn force(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_force(&(*potential_distance/self.contour_length), &(potential_stiffness*self.contour_length.powf(2.0)/BOLTZMANN_CONSTANT/temperature))*BOLTZMANN_CONSTANT*temperature/self.link_length
    }
    fn nondimensional_force(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        let integrand_numerator = |nondimensional_end_to_end_length_per_link: f64|
        {
            let mut sum: f64 = 0.0;
            let n = self.number_of_links as u128;
            let p = self.number_of_links_f64 - 2.0;
            let m = -nondimensional_end_to_end_length_per_link*0.5 + 0.5;
            let k = (self.number_of_links_f64*m).ceil() as u128;
            for s in 0..k
            {
                sum += (-1.0_f64).powf(s as f64)*((factorial(n)/factorial(s)/factorial(n - s)) as f64)*(m - (s as f64)/self.number_of_links_f64).powf(p);
            }
            0.5*nondimensional_end_to_end_length_per_link*(n.pow(n as u32) as f64)/(factorial(n - 2) as f64)*sum*((nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link)*(-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link).powf(2.0)).exp() - nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link)*(-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link).powf(2.0)).exp())/(2.0*nondimensional_potential_stiffness*nondimensional_potential_distance*nondimensional_end_to_end_length_per_link) + ((-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link).powf(2.0)).exp() - (-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link).powf(2.0)).exp())/(2.0*nondimensional_potential_stiffness*nondimensional_potential_distance.powf(2.0)*nondimensional_end_to_end_length_per_link))
        };
        let integrand_denominator = |nondimensional_end_to_end_length_per_link: f64|
        {
            let mut sum: f64 = 0.0;
            let n = self.number_of_links as u128;
            let p = self.number_of_links_f64 - 2.0;
            let m = -nondimensional_end_to_end_length_per_link*0.5 + 0.5;
            let k = (self.number_of_links_f64*m).ceil() as u128;
            for s in 0..k
            {
                sum += (-1.0_f64).powf(s as f64)*((factorial(n)/factorial(s)/factorial(n - s)) as f64)*(m - (s as f64)/self.number_of_links_f64).powf(p);
            }
            0.5*nondimensional_end_to_end_length_per_link*(n.pow(n as u32) as f64)/(factorial(n - 2) as f64)*sum*((-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link).powf(2.0)).exp() - (-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link).powf(2.0)).exp())/(2.0*nondimensional_potential_stiffness*nondimensional_potential_distance*nondimensional_end_to_end_length_per_link)
        };
        integrate(integrand_numerator, ZERO, ONE, POINTS)/integrate(integrand_denominator, ZERO, ONE, POINTS)/self.number_of_links_f64
    
    }
    fn helmholtz_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(&(potential_distance/self.contour_length), &(potential_stiffness*self.contour_length.powf(2.0)/BOLTZMANN_CONSTANT/temperature), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    fn helmholtz_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy_per_link(&(*potential_distance/self.contour_length), &(potential_stiffness*self.contour_length.powf(2.0)/BOLTZMANN_CONSTANT/temperature), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    fn relative_helmholtz_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy(&(*potential_distance/self.contour_length), &(potential_stiffness*self.contour_length.powf(2.0)/BOLTZMANN_CONSTANT/temperature))*BOLTZMANN_CONSTANT*temperature
    }
    fn relative_helmholtz_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(&(*potential_distance/self.contour_length), &(potential_stiffness*self.contour_length.powf(2.0)/BOLTZMANN_CONSTANT/temperature))*BOLTZMANN_CONSTANT*temperature
    }
    fn nondimensional_helmholtz_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let integrand = |nondimensional_end_to_end_length_per_link: f64|
        {
            let mut sum: f64 = 0.0;
            let n = self.number_of_links as u128;
            let p = self.number_of_links_f64 - 2.0;
            let m = -nondimensional_end_to_end_length_per_link*0.5 + 0.5;
            let k = (self.number_of_links_f64*m).ceil() as u128;
            for s in 0..k
            {
                sum += (-1.0_f64).powf(s as f64)*((factorial(n)/factorial(s)/factorial(n - s)) as f64)*(m - (s as f64)/self.number_of_links_f64).powf(p);
            }
            0.5*nondimensional_end_to_end_length_per_link*(n.pow(n as u32) as f64)/(factorial(n - 2) as f64)*sum*((-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link).powf(2.0)).exp() - (-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link).powf(2.0)).exp())/(2.0*nondimensional_potential_stiffness*nondimensional_potential_distance*nondimensional_end_to_end_length_per_link)
        };
        let nondimensional_configurational_partition_function = integrate(integrand, ZERO, ONE, POINTS);
        -ln(&nondimensional_configurational_partition_function) - self.number_of_links_f64*ln(&(8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)))
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
    fn gibbs_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.helmholtz_free_energy(potential_distance, potential_stiffness, temperature) - 0.5*potential_stiffness*potential_distance.powf(2.0)
    }
    fn gibbs_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.helmholtz_free_energy_per_link(potential_distance, potential_stiffness, temperature) - 0.5*potential_stiffness*potential_distance.powf(2.0)/self.number_of_links_f64
    }
    fn relative_gibbs_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.relative_helmholtz_free_energy(potential_distance, potential_stiffness, temperature) - 0.5*potential_stiffness*potential_distance.powf(2.0)
    }
    fn relative_gibbs_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.relative_helmholtz_free_energy_per_link(potential_distance, potential_stiffness, temperature) - 0.5*potential_stiffness*potential_distance.powf(2.0)/self.number_of_links_f64
    }
    fn nondimensional_gibbs_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness, temperature) - 0.5*nondimensional_potential_stiffness*nondimensional_potential_distance.powf(2.0)
    }
    fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy_per_link(nondimensional_potential_distance, nondimensional_potential_stiffness, temperature) - 0.5*nondimensional_potential_stiffness*nondimensional_potential_distance.powf(2.0)/self.number_of_links_f64
    }
    fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness) - 0.5*nondimensional_potential_stiffness*nondimensional_potential_distance.powf(2.0)
    }
    fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_potential_distance, nondimensional_potential_stiffness) - 0.5*nondimensional_potential_stiffness*nondimensional_potential_distance.powf(2.0)/self.number_of_links_f64
    }
}