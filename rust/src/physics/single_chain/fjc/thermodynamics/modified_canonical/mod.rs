use std::f64::consts::PI;
use crate::math::
{
    ln,
    factorial,
    sinhc,
    integrate
};
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::ModifiedCanonical;

pub mod test;

static ZERO: f64 = 1e-8;

pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u16,
    pub number_of_links_f64: f64,
    pub contour_length: f64,
    pub potential_stiffness: f64
}

impl ModifiedCanonical for FJC
{
    fn init(number_of_links: u16, link_length: f64, hinge_mass: f64, potential_stiffness: f64) -> FJC
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length,
            potential_stiffness: potential_stiffness
        }
    }
    fn helmholtz_free_energy(&self, potential_distance: &f64, temperature: f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(&(*potential_distance/self.contour_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    fn helmholtz_free_energy_per_link(&self, potential_distance: &f64, temperature: f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy_per_link(&(*potential_distance/self.contour_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    fn relative_helmholtz_free_energy(&self, potential_distance: &f64, temperature: f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(&(*potential_distance/self.contour_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    fn relative_helmholtz_free_energy_per_link(&self, potential_distance: &f64, temperature: f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(&(*potential_distance/self.contour_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    fn nondimensional_helmholtz_free_energy(&self, nondimensional_potential_distance: &f64, temperature: f64) -> f64
    {
        let nondimensional_potential_stiffness = BOLTZMANN_CONSTANT*temperature*self.potential_stiffness*self.contour_length.powf(2.0);
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
            0.5*nondimensional_end_to_end_length_per_link*(n.pow(n as u32) as f64)/(factorial(n - 2) as f64)*sum*0.5*((-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link).powf(2.0)).exp() - (-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link).powf(2.0)).exp())/(nondimensional_potential_stiffness*nondimensional_potential_distance*nondimensional_end_to_end_length_per_link)
        };
        let nondimensional_configurational_partition_function = integrate(integrand, 0.0, 1.0, 10000);
        -ln(&nondimensional_configurational_partition_function) - self.number_of_links_f64*ln(&(8.0*PI.powf(2.0)*self.hinge_mass*self.link_length.powf(2.0)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powf(2.0)))
    }
    fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_potential_distance: &f64, temperature: f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(nondimensional_potential_distance, temperature)/self.number_of_links_f64
    }
    fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_potential_distance: &f64, temperature: f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(nondimensional_potential_distance, temperature) - self.nondimensional_helmholtz_free_energy(&ZERO, temperature)
    }
    fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_potential_distance: &f64, temperature: f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy(nondimensional_potential_distance, temperature)/self.number_of_links_f64
    }
}