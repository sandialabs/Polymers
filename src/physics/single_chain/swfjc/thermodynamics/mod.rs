mod test;
pub mod isotensional;
pub struct SWFJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub well_width: f64,
    pub isotensional: self::isotensional::SWFJC,
}
impl SWFJC
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, well_width: f64) -> Self
    {
        SWFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            well_width,
            isotensional: self::isotensional::SWFJC::init(number_of_links, link_length, hinge_mass, well_width),
        }
    }
}
