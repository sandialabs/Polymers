use crate::physics::single_chain::
{
    Isometric,
    // Isotensional
};

pub mod test;
pub mod isometric;
pub mod isotensional;

pub struct Thermodynamics
{
    pub link_length: f64,
    pub number_of_links: u16,
    pub isometric: isometric::FJC,
    pub isotensional: isotensional::Isotensional
}

impl Thermodynamics
{
    pub fn init(number_of_links: u16, link_length: f64) -> Thermodynamics
    {
        Thermodynamics
        {
            link_length: link_length,
            number_of_links: number_of_links,
            isometric: isometric::FJC::init(number_of_links, link_length),
            isotensional: isotensional::Isotensional::init(number_of_links, link_length)
        }
    }
}
