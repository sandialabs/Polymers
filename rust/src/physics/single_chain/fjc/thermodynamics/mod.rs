use crate::physics::single_chain::fjc::thermodynamics::isometric::Isometric;
use crate::physics::single_chain::fjc::thermodynamics::isotensional::Isotensional;

pub mod test;
pub mod isometric;
pub mod isotensional;

pub struct Thermodynamics
{
    pub number_of_links: u16,
    pub isometric: Isometric,
    pub isotensional: Isotensional
}

impl Thermodynamics
{
    pub fn init(number_of_links: u16) -> Thermodynamics
    {
        Thermodynamics
        {
            number_of_links: number_of_links,
            isometric: Isometric::init(number_of_links),
            isotensional: Isotensional::init(number_of_links)
        }
    }
}