pub mod test;
pub mod mechanics;
pub mod thermodynamics;

// move to single_chain?
//pub static ONE: f64 = 1.0;
//pub static ZERO: f64 = 1e-6;
//pub static POINTS: u128 = 100;

pub struct SWFJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub well_width: f64,
    pub mechanics: mechanics::SWFJC,
    pub thermodynamics: thermodynamics::SWFJC
}
impl SWFJC
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, well_width: f64) -> SWFJC
    {
        SWFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            well_width,
            mechanics: mechanics::SWFJC::init(number_of_links, link_length, hinge_mass, well_width),
            thermodynamics: thermodynamics::SWFJC::init(number_of_links, link_length, hinge_mass, well_width),
        }
    }
}
