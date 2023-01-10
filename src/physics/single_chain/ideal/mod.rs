mod test;
pub mod thermodynamics;
pub struct Ideal
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub thermodynamics: thermodynamics::Ideal
}
impl Ideal
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
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
