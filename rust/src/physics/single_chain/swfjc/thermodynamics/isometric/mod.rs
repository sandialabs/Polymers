pub mod test;
pub mod legendre;
use crate::physics::single_chain::swfjc::thermodynamics::isometric::legendre::SWFJC as SWFJCLegendre;
use crate::physics::single_chain::swfjc::thermodynamics::
{
    Isometric,
    IsometricLegendre
};
pub struct SWFJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub number_of_links_f64: f64,
    pub well_width: f64,
    pub contour_length: f64,
    pub nondimensional_well_parameter: f64,
    pub legendre: SWFJCLegendre
}
impl Isometric for SWFJC
{
    fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, well_width: f64) -> SWFJC
    {
        SWFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            well_width,
            number_of_links_f64: number_of_links as f64,
            contour_length: (number_of_links as f64)*link_length,
            nondimensional_well_parameter: 1.0 + well_width/link_length,
            legendre: SWFJCLegendre::init(number_of_links, link_length, hinge_mass, well_width)
        }
    }
}
