#[cfg(feature = "python")]
pub mod py;

mod test;

/// ????????????????????????????????
pub mod weak_potential;

/// ????????????????????????????????
pub mod strong_potential;

/// The structure of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach.
pub struct FJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// ????????????????????????????????
    pub weak_potential: weak_potential::FJC,

    /// ????????????????????????????????
    pub strong_potential: strong_potential::FJC
}

/// The implemented functionality of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach.
impl FJC
{
    /// Initializes and returns an instance of the thermodynamics of the FJC model in the modified canonical ensemble approximated using an asymptotic approach.
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
