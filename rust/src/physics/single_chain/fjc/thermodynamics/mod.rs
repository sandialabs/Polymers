use crate::physics::single_chain::
{
    Isometric,
    Isotensional
};

pub mod test;
pub mod isometric;
pub mod isotensional;

pub struct Thermodynamics
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u16,
    pub isometric: isometric::FJC,
    pub isotensional: isotensional::FJC
}

impl Thermodynamics
{
    pub fn init(number_of_links: u16, link_length: f64, hinge_mass: f64) -> Thermodynamics
    {
        Thermodynamics
        {
            hinge_mass,
            link_length,
            number_of_links,
            isometric: isometric::FJC::init(number_of_links, link_length, hinge_mass),
            isotensional: isotensional::FJC::init(number_of_links, link_length, hinge_mass)
        }
    }
}
