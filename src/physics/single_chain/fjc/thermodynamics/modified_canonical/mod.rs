mod test;
pub mod asymptotic;
use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::
{
    ONE,
    ZERO,
    POINTS
};
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    number_of_links_f64: f64,
    contour_length: f64,
    pub asymptotic: asymptotic::FJC,
}
impl FJC
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length,
            asymptotic: asymptotic::FJC::init(number_of_links, link_length, hinge_mass)
        }
    }
    pub fn end_to_end_length(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        potential_distance - self.force(potential_distance, potential_stiffness, temperature)/potential_stiffness
    }
    pub fn end_to_end_length_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        (potential_distance - self.force(potential_distance, potential_stiffness, temperature)/potential_stiffness)/self.number_of_links_f64
    }
    pub fn nondimensional_end_to_end_length(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        self.number_of_links_f64*nondimensional_potential_distance - self.number_of_links_f64.powi(2)*self.nondimensional_force(nondimensional_potential_distance, nondimensional_potential_stiffness)/nondimensional_potential_stiffness
    }
    pub fn nondimensional_end_to_end_length_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        nondimensional_potential_distance - self.number_of_links_f64*self.nondimensional_force(nondimensional_potential_distance, nondimensional_potential_stiffness)/nondimensional_potential_stiffness
    }
    pub fn force(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_force(&(*potential_distance/self.contour_length), &(potential_stiffness*self.contour_length.powi(2)/BOLTZMANN_CONSTANT/temperature))*BOLTZMANN_CONSTANT*temperature/self.link_length
    }
    pub fn nondimensional_force(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        let n = self.number_of_links as u128;
        let p: i32 = (self.number_of_links - 2).into();
        let integrand_numerator = |nondimensional_end_to_end_length_per_link: f64|
        {
            let m = -nondimensional_end_to_end_length_per_link*0.5 + 0.5;
            let k = (self.number_of_links_f64*m).ceil() as u128;
            let sum: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powi(p)).sum();
            0.5*nondimensional_end_to_end_length_per_link*(n.pow(n as u32) as f64)/((1..=n-2).product::<u128>() as f64)*sum*((nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link)*(-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link).powi(2)).exp() - nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link)*(-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link).powi(2)).exp())/(2.0*nondimensional_potential_stiffness*nondimensional_potential_distance*nondimensional_end_to_end_length_per_link) + ((-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link).powi(2)).exp() - (-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link).powi(2)).exp())/(2.0*nondimensional_potential_stiffness*nondimensional_potential_distance.powi(2)*nondimensional_end_to_end_length_per_link))
        };
        let integrand_denominator = |nondimensional_end_to_end_length_per_link: f64|
        {
            let m = -nondimensional_end_to_end_length_per_link*0.5 + 0.5;
            let k = (self.number_of_links_f64*m).ceil() as u128;
            let sum: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powi(p)).sum();
            0.5*nondimensional_end_to_end_length_per_link*(n.pow(n as u32) as f64)/((1..=n-2).product::<u128>() as f64)*sum*((-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link).powi(2)).exp() - (-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link).powi(2)).exp())/(2.0*nondimensional_potential_stiffness*nondimensional_potential_distance*nondimensional_end_to_end_length_per_link)
        };
        let dx = (ONE - ZERO)/(POINTS as f64);
        (0..=POINTS-1).collect::<Vec::<u128>>().iter().map(|index| integrand_numerator(ZERO + (0.5 + *index as f64)*dx)).sum::<f64>()/(0..=POINTS-1).collect::<Vec::<u128>>().iter().map(|index| integrand_denominator(ZERO + (0.5 + *index as f64)*dx)).sum::<f64>()/self.number_of_links_f64    
    }
    pub fn helmholtz_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(&(potential_distance/self.contour_length), &(potential_stiffness*self.contour_length.powi(2)/BOLTZMANN_CONSTANT/temperature), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    pub fn helmholtz_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy_per_link(&(*potential_distance/self.contour_length), &(potential_stiffness*self.contour_length.powi(2)/BOLTZMANN_CONSTANT/temperature), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    pub fn relative_helmholtz_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy(&(*potential_distance/self.contour_length), &(potential_stiffness*self.contour_length.powi(2)/BOLTZMANN_CONSTANT/temperature))*BOLTZMANN_CONSTANT*temperature
    }
    pub fn relative_helmholtz_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(&(*potential_distance/self.contour_length), &(potential_stiffness*self.contour_length.powi(2)/BOLTZMANN_CONSTANT/temperature))*BOLTZMANN_CONSTANT*temperature
    }
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let integrand = |nondimensional_end_to_end_length_per_link: f64|
        {
            let n = self.number_of_links as u128;
            let p: i32 = (self.number_of_links - 2).into();
            let m = -nondimensional_end_to_end_length_per_link*0.5 + 0.5;
            let k = (self.number_of_links_f64*m).ceil() as u128;
            let sum: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powi(p)).sum();
            0.5*nondimensional_end_to_end_length_per_link*(n.pow(n as u32) as f64)/((1..=n-2).product::<u128>() as f64)*sum*((-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance - nondimensional_end_to_end_length_per_link).powi(2)).exp() - (-0.5*nondimensional_potential_stiffness*(nondimensional_potential_distance + nondimensional_end_to_end_length_per_link).powi(2)).exp())/(2.0*nondimensional_potential_stiffness*nondimensional_potential_distance*nondimensional_end_to_end_length_per_link)
        };
        let dx = (ONE - ZERO)/(POINTS as f64);
        let nondimensional_configurational_partition_function = (0..=POINTS-1).collect::<Vec::<u128>>().iter().map(|index| integrand(ZERO + (0.5 + *index as f64)*dx)).sum::<f64>()*dx;
        -nondimensional_configurational_partition_function.ln() - (self.number_of_links_f64 - 1.0)*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness, temperature)/self.number_of_links_f64
    }
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness, &300.0) - self.nondimensional_helmholtz_free_energy(&ZERO, nondimensional_potential_stiffness, &300.0)
    }
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness)/self.number_of_links_f64
    }
    pub fn gibbs_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.helmholtz_free_energy(potential_distance, potential_stiffness, temperature) - 0.5*potential_stiffness*potential_distance.powi(2)
    }
    pub fn gibbs_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.helmholtz_free_energy_per_link(potential_distance, potential_stiffness, temperature) - 0.5*potential_stiffness*potential_distance.powi(2)/self.number_of_links_f64
    }
    pub fn relative_gibbs_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.relative_helmholtz_free_energy(potential_distance, potential_stiffness, temperature) - 0.5*potential_stiffness*potential_distance.powi(2)
    }
    pub fn relative_gibbs_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.relative_helmholtz_free_energy_per_link(potential_distance, potential_stiffness, temperature) - 0.5*potential_stiffness*potential_distance.powi(2)/self.number_of_links_f64
    }
    pub fn nondimensional_gibbs_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness, temperature) - 0.5*nondimensional_potential_stiffness*nondimensional_potential_distance.powi(2)
    }
    pub fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_helmholtz_free_energy_per_link(nondimensional_potential_distance, nondimensional_potential_stiffness, temperature) - 0.5*nondimensional_potential_stiffness*nondimensional_potential_distance.powi(2)/self.number_of_links_f64
    }
    pub fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy(nondimensional_potential_distance, nondimensional_potential_stiffness) - 0.5*nondimensional_potential_stiffness*nondimensional_potential_distance.powi(2)
    }
    pub fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_potential_distance, nondimensional_potential_stiffness) - 0.5*nondimensional_potential_stiffness*nondimensional_potential_distance.powi(2)/self.number_of_links_f64
    }
}
