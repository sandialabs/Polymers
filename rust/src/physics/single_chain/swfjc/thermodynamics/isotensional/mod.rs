pub mod test;
use crate::physics::single_chain::swfjc::thermodynamics::isotensional::legendre::SWFJC as SWFJCLegendre;
use crate::physics::single_chain::swfjc::thermodynamics::
{
    Isotensional,
    IsotensionalLegendre
};
pub struct SWFJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub well_width: f64,
    pub legendre: SWFJCLegendre
}
impl Isotensional for SWFJC
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, well_width: f64) -> SWFJC
    {
        SWFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            well_width,
            legendre: SWFJCLegendre::init(number_of_links, link_length, hinge_mass, well_width)
        }
    }
}
