pub mod test;
pub mod strong_potential;
pub mod weak_potential;
use crate::physics::single_chain::fjc::thermodynamics::{
    ModifiedCanonicalAsymptotic,
    ModifiedCanonicalAsymptoticStrongPotential,
    ModifiedCanonicalAsymptoticWeakPotential
};
use crate::physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::strong_potential::FJC as FJCAsymptoticStrongPotential;
use crate::physics::single_chain::fjc::thermodynamics::modified_canonical::asymptotic::weak_potential::FJC as FJCAsymptoticWeakPotential;
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub number_of_links_f64: f64,
    pub contour_length: f64,
    pub strong_potential: FJCAsymptoticStrongPotential,
    pub weak_potential: FJCAsymptoticWeakPotential
}
impl ModifiedCanonicalAsymptotic for FJC
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
            strong_potential: FJCAsymptoticStrongPotential::init(number_of_links, link_length, hinge_mass),
            weak_potential: FJCAsymptoticWeakPotential::init(number_of_links, link_length, hinge_mass)
        }
    }
}