pub mod test;
pub mod thermodynamics;
use self::
{
    thermodynamics::
    {
        Isometric,
        Isotensional,
        ModifiedCanonical,
        isometric::Legendre as IsometricLegendre,
        isotensional::Legendre as IsotensionalLegendre,
        modified_canonical::
        {
            asymptotic::
            {
                WeakPotential as ModifiedCanonicalAsymptoticWeakPotential,
                StrongPotential as ModifiedCanonicalAsymptoticStrongPotential
            }
        }
    }
};
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub thermodynamics: thermodynamics::FJC
}
impl FJC
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> FJC
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            thermodynamics: thermodynamics::FJC::init(number_of_links, link_length, hinge_mass),
        }
    }
}
