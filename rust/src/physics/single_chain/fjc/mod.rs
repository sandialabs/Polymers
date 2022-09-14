pub mod test;
pub mod mechanics;
pub mod thermodynamics;

pub struct FJC
{
    pub link_length: f64,
    pub number_of_links: u16,
    pub mechanics: mechanics::Mechanics,
    pub thermodynamics: thermodynamics::Thermodynamics
}

impl FJC
{
    pub fn init(number_of_links: u16, link_length: f64) -> FJC
    {
        FJC
        {
            link_length: link_length,
            number_of_links: number_of_links,
            mechanics: mechanics::Mechanics::init(number_of_links, link_length),
            thermodynamics: thermodynamics::Thermodynamics::init(number_of_links, link_length),
        }
    }
}
