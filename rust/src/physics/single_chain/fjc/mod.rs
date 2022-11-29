pub mod test;
pub mod mechanics;
pub mod thermodynamics;
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub mechanics: mechanics::Mechanics,
    pub thermodynamics: thermodynamics::Thermodynamics
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
            mechanics: mechanics::Mechanics::init(number_of_links, link_length, hinge_mass),
            thermodynamics: thermodynamics::Thermodynamics::init(number_of_links, link_length, hinge_mass),
        }
    }
}