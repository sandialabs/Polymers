mod test;
pub mod weak_potential;
pub mod strong_potential;
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
    pub weak_potential: weak_potential::FJC,
    pub strong_potential: strong_potential::FJC
}
impl FJC
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            weak_potential: weak_potential::FJC::init(number_of_links, link_length, hinge_mass),
            strong_potential: strong_potential::FJC::init(number_of_links, link_length, hinge_mass)
        }
    }
}
