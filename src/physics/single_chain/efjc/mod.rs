mod test;
pub mod thermodynamics;
pub struct EFJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub link_stiffness: f64,
    pub thermodynamics: self::thermodynamics::EFJC
}
impl EFJC
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        EFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            thermodynamics: self::thermodynamics::EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness),
        }
    }
}
