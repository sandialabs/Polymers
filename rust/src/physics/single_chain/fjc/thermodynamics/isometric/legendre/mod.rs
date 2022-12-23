pub mod test;
use std::f64::consts::PI;
use crate::math::integrate;
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::fjc::
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
    pub number_of_links_f64: f64,
    pub contour_length: f64,
    normalization_nondimensional_equilibrium_distribution: f64
}
use super::Legendre;
impl Legendre for FJC
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> FJC
    {
        let temporary_model = FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length,
            normalization_nondimensional_equilibrium_distribution: 1.0
        };
        let dx = (ONE - ZERO)/(POINTS as f64);
        let normalization = (0..=POINTS-1).collect::<Vec::<u128>>().iter().map(|index| temporary_model.nondimensional_equilibrium_radial_distribution(&(ZERO + (0.5 + *index as f64)*dx))).sum::<f64>()*dx;
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length,
            normalization_nondimensional_equilibrium_distribution: normalization
        }
    }
    fn force(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        let nondimensional_end_to_end_length_per_link = end_to_end_length/self.contour_length;
        BOLTZMANN_CONSTANT*temperature/self.link_length*(2.14234*nondimensional_end_to_end_length_per_link.powi(3) - 4.22785*nondimensional_end_to_end_length_per_link.powi(2) + 3.0*nondimensional_end_to_end_length_per_link)/(1.0 - nondimensional_end_to_end_length_per_link)/(0.71716*nondimensional_end_to_end_length_per_link.powi(3) - 0.41103*nondimensional_end_to_end_length_per_link.powi(2) - 0.39165*nondimensional_end_to_end_length_per_link + 1.0)
    }
    fn nondimensional_force(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        (2.14234*nondimensional_end_to_end_length_per_link.powi(3) - 4.22785*nondimensional_end_to_end_length_per_link.powi(2) + 3.0*nondimensional_end_to_end_length_per_link)/(1.0 - nondimensional_end_to_end_length_per_link)/(0.71716*nondimensional_end_to_end_length_per_link.powi(3) - 0.41103*nondimensional_end_to_end_length_per_link.powi(2) - 0.39165*nondimensional_end_to_end_length_per_link + 1.0)
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
        self.nondimensional_relative_helmholtz_free_energy(nondimensional_end_to_end_length_per_link) - (self.number_of_links_f64 - 1.0)*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_end_to_end_length_per_link) - (self.number_of_links_f64 - 1.0)/self.number_of_links_f64*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        self.nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_end_to_end_length_per_link)*self.number_of_links_f64
    }
    fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        let nondimensional_force = self.nondimensional_force(nondimensional_end_to_end_length_per_link);
        nondimensional_force**nondimensional_end_to_end_length_per_link - (nondimensional_force.sinh()/nondimensional_force).ln()
    }
    fn equilibrium_distribution(&self, end_to_end_length: &f64) -> f64
    {
        self.nondimensional_equilibrium_distribution(&(end_to_end_length/self.contour_length))/self.contour_length.powi(3)
    }
    fn nondimensional_equilibrium_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        let nondimensional_force = self.nondimensional_force(nondimensional_end_to_end_length_per_link);
        (nondimensional_force.sinh()/nondimensional_force*(-nondimensional_force**nondimensional_end_to_end_length_per_link).exp()).powi(self.number_of_links as i32)/self.normalization_nondimensional_equilibrium_distribution
    }
    fn equilibrium_radial_distribution(&self, end_to_end_length: &f64) -> f64
    {
        self.nondimensional_equilibrium_radial_distribution(&(end_to_end_length/self.contour_length))/self.contour_length
    }
    fn nondimensional_equilibrium_radial_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        4.0*PI*nondimensional_end_to_end_length_per_link.powi(2)*self.nondimensional_equilibrium_distribution(nondimensional_end_to_end_length_per_link)
    }
    fn gibbs_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_gibbs_free_energy(&(end_to_end_length/self.contour_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    fn gibbs_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_gibbs_free_energy_per_link(&(end_to_end_length/self.contour_length), temperature)*BOLTZMANN_CONSTANT*temperature
    }
    fn relative_gibbs_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        self.relative_gibbs_free_energy_per_link(end_to_end_length, temperature)*self.number_of_links_f64
    }
    fn relative_gibbs_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_relative_gibbs_free_energy_per_link(&(end_to_end_length/self.contour_length))*BOLTZMANN_CONSTANT*temperature
    }
    fn nondimensional_gibbs_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        let n = self.number_of_links as u128;
        let p: i32 = (self.number_of_links - 2).into();
        let m = -*nondimensional_end_to_end_length_per_link*0.5 + 0.5;
        let k = (self.number_of_links_f64*m).ceil() as u128;
        let sum: f64 = (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/self.number_of_links_f64).powi(p)).sum();
        -(nondimensional_end_to_end_length_per_link*self.number_of_links_f64)*self.nondimensional_force(nondimensional_end_to_end_length_per_link) - (0.125/PI/nondimensional_end_to_end_length_per_link*(n.pow(n as u32) as f64)/((1..=n-2).product::<u128>() as f64)*sum/self.contour_length.powi(3)).ln() - (self.number_of_links_f64 - 1.0)*(8.0*PI.powi(2)*self.hinge_mass*self.link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
    }
    fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        self.nondimensional_gibbs_free_energy(nondimensional_end_to_end_length_per_link, temperature)/self.number_of_links_f64
    }
    fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        self.nondimensional_gibbs_free_energy(nondimensional_end_to_end_length_per_link, &300.0) - self.nondimensional_gibbs_free_energy(&ZERO, &300.0)
    }
    fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        self.nondimensional_relative_gibbs_free_energy(nondimensional_end_to_end_length_per_link)/self.number_of_links_f64
    }
}
