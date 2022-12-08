pub mod test;
pub mod mechanics;
pub mod thermodynamics;
pub static ONE: f64 = 1.0;
pub static ZERO: f64 = 1e-6;
pub static POINTS: u128 = 100;
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub mechanics: mechanics::FJC,
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
            mechanics: mechanics::FJC::init(number_of_links, link_length, hinge_mass),
            thermodynamics: thermodynamics::FJC::init(number_of_links, link_length, hinge_mass),
        }
    }
}