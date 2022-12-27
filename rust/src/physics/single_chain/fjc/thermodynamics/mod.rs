pub mod test;
pub mod isometric;
pub mod isotensional;
pub mod modified_canonical;
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub isometric: isometric::FJC,
    pub isotensional: isotensional::FJC,
    pub modified_canonical: modified_canonical::FJC
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
            isometric: isometric::FJC::init(number_of_links, link_length, hinge_mass),
            isotensional: isotensional::FJC::init(number_of_links, link_length, hinge_mass),
            modified_canonical: modified_canonical::FJC::init(number_of_links, link_length, hinge_mass)
        }
    }
}
