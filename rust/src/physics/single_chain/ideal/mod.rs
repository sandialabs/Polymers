pub mod test;
pub mod thermodynamics;
pub static ONE: f64 = 1.0;
pub static ZERO: f64 = 1e-6;
pub static POINTS: u128 = 100;
use self::
{
    thermodynamics::
    {
        Isometric,
        Isotensional
    }
};
pub struct Ideal
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub thermodynamics: thermodynamics::Ideal
}
impl Ideal
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Ideal
    {
        Ideal
        {
            hinge_mass,
            link_length,
            number_of_links,
            thermodynamics: thermodynamics::Ideal::init(number_of_links, link_length, hinge_mass),
        }
    }
}
