#[cfg(feature = "python")]
pub mod py;

mod test;

/// The freely-jointed chain (FJC) model thermodynamics.
pub mod thermodynamics;

/// The structure of the FJC model.
pub struct FJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The thermodynamic functions of the model.
    pub thermodynamics: thermodynamics::FJC
}

/// The implemented functionality of the FJC model.
impl FJC
{
    /// Initializes and returns an instance of the FJC model.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            thermodynamics: thermodynamics::FJC::init(number_of_links, link_length, hinge_mass),
        }
    }
}
