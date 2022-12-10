pub mod test;
pub struct SWFJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
}
impl SWFJC
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> SWFJC
    {
        SWFJC
        {
            hinge_mass,
            link_length,
            number_of_links
        }
    }
}
