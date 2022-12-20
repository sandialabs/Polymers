pub mod test;
pub mod isotensional;
use self::isotensional::{
    Legendre as IsotensionalLegendre,
    asymptotic::{
        Alternative as IsotensionalAsymptoticAlternative,
        Reduced as IsotensionalAsymptoticReduced,
        alternative::Legendre as IsotensionalAsymptoticAlternativeLegendre,
        reduced::Legendre as IsotensionalAsymptoticReducedLegendre
    }
};
pub struct EFJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub link_stiffness: f64,
    pub isotensional: self::isotensional::EFJC,
}
impl EFJC
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> EFJC
    {
        EFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            isotensional: self::isotensional::EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness),
        }
    }
}
pub trait Isotensional
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
