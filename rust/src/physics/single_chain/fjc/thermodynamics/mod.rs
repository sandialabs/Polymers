use crate::physics::single_chain::fjc::thermodynamics::isometric::Isometric;
use crate::physics::single_chain::fjc::thermodynamics::isotensional::Isotensional;

pub mod test;
pub mod isometric;
pub mod isotensional;

pub struct Thermodynamics
{
    pub number_of_links: u16,
    pub link_length: f64,
    pub isometric: Isometric,
    pub isotensional: Isotensional
}

impl Thermodynamics
{
    pub fn init(number_of_links: u16, link_length: f64) -> Thermodynamics
    {
        Thermodynamics
        {
            number_of_links: number_of_links,
            link_length: link_length,
            isometric: Isometric::init(number_of_links, link_length),
            isotensional: Isotensional::init(number_of_links, link_length)
        }
    }
}
