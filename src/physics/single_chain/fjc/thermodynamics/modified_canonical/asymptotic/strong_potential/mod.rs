mod test;
use std::f64::consts::PI;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::ZERO;
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    number_of_links_f64: f64,
    contour_length: f64
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
            contour_length: (number_of_links as f64)*link_length
        }
    }
    pub fn force(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        BOLTZMANN_CONSTANT*temperature/self.link_length*self.nondimensional_force(&(potential_distance/self.contour_length), &(potential_stiffness*(self.contour_length).powi(2)/BOLTZMANN_CONSTANT/temperature))
    }
    pub fn nondimensional_force(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64) -> f64
    {
        let n = self.number_of_links as u128;
        let p: i32 = (self.number_of_links - 2).into();
        let m = -*nondimensional_potential_distance*0.5 + 0.5;
        let k = (self.number_of_links_f64*m).ceil() as u128;
        let sum_0: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powi(p)).sum();
        let sum_1: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powi(p - 1)).sum();
        let sum_2: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powi(p - 2)).sum();
        let sum_3: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powi(p - 3)).sum();
        (1.0/nondimensional_potential_distance + (0.5*self.number_of_links_f64 - 1.0)*sum_1/sum_0)/self.number_of_links_f64 + 0.5/nondimensional_potential_stiffness/self.number_of_links_f64*((0.5*self.number_of_links_f64 - 1.0)*((0.5*self.number_of_links_f64 - 1.0)*sum_1/sum_0*((self.number_of_links_f64 - 2.0)*(sum_1/sum_0).powi(2) - (self.number_of_links_f64 - 3.0)*sum_2/sum_0) - (0.5*self.number_of_links_f64 - 1.5)*((0.5*self.number_of_links_f64 - 1.0)*sum_1*sum_2/sum_0.powi(2) - (0.5*self.number_of_links_f64 - 2.0)*sum_3/sum_0)) + 2.0*nondimensional_potential_distance.powi(-3) - 2.0*((0.5*self.number_of_links_f64 - 1.0)*sum_1/sum_0 + nondimensional_potential_distance.powi(-1))*((0.5*self.number_of_links_f64 - 1.0)*((0.5*self.number_of_links_f64 - 1.0)*(sum_1/sum_0).powi(2) - (0.5*self.number_of_links_f64 - 1.5)*sum_2/sum_0) - nondimensional_potential_distance.powi(-2)))
    }
    pub fn helmholtz_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        BOLTZMANN_CONSTANT*temperature*self.nondimensional_helmholtz_free_energy(&(potential_distance/self.contour_length), &(potential_stiffness*(self.contour_length).powi(2)/BOLTZMANN_CONSTANT/temperature), temperature)
    }
    pub fn helmholtz_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.helmholtz_free_energy(potential_distance, potential_stiffness, temperature)/self.number_of_links_f64
    }
    pub fn relative_helmholtz_free_energy(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.helmholtz_free_energy(potential_distance, potential_stiffness, temperature) - self.helmholtz_free_energy(&(ZERO*self.number_of_links_f64*self.link_length), potential_stiffness, temperature)
    }
    pub fn relative_helmholtz_free_energy_per_link(&self, potential_distance: &f64, potential_stiffness: &f64, temperature: &f64) -> f64
    {
        self.relative_helmholtz_free_energy(potential_distance, potential_stiffness, temperature)/self.number_of_links_f64
    }
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_potential_distance: &f64, nondimensional_potential_stiffness: &f64, temperature: &f64) -> f64
    {
        let n = self.number_of_links as u128;
        let p: i32 = (self.number_of_links - 2).into();
        let m = -*nondimensional_potential_distance*0.5 + 0.5;
        let k = (self.number_of_links_f64*m).ceil() as u128;
        let sum_0: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powi(p)).sum();
        let sum_1: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powi(p - 1)).sum();
        let sum_2: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powi(p - 2)).sum();
        -(0.125/PI/nondimensional_potential_distance*(n.pow(n as u32) as f64)/((1..=n-2).product::<u128>() as f64)*sum_0/self.contour_length.powi(3)).ln() - (self.number_of_links_f64 - 1.0)*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln() - 1.5*(2.0*PI/nondimensional_potential_stiffness).ln() - 3.0*(self.contour_length).ln() + 0.5/nondimensional_potential_stiffness*((0.5*self.number_of_links_f64 - 1.0)*((0.5*self.number_of_links_f64 - 1.0)*(sum_1/sum_0).powi(2) - (0.5*self.number_of_links_f64 - 1.5)*sum_2/sum_0) - nondimensional_potential_distance.powi(-2) - ((0.5*self.number_of_links_f64 - 1.0)*sum_1/sum_0 + nondimensional_potential_distance.powi(-1)).powi(2))
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
}
