mod test;
pub mod isometric;
pub mod isotensional;
pub struct Ideal
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub isometric: isometric::Ideal,
    pub isotensional: isotensional::Ideal
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
            isometric: isometric::Ideal::init(number_of_links, link_length, hinge_mass),
            isotensional: isotensional::Ideal::init(number_of_links, link_length, hinge_mass)
        }
    }
}
