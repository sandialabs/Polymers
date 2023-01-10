#[cfg(feature = "python")]
pub mod py;

mod test;

/// The freely-jointed chain (FJC) model thermodynamics in the isometric ensemble.
pub mod isometric;

/// The freely-jointed chain (FJC) model thermodynamics in the isotensional ensemble.
pub mod isotensional;

/// The freely-jointed chain (FJC) model thermodynamics in the modified canonical ensemble.
pub mod modified_canonical;

/// The structure of the thermodynamics of the FJC model.
pub struct FJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The thermodynamic functions of the model in the isometric ensemble.
    pub isometric: isometric::FJC,

    /// The thermodynamic functions of the model in the isotensional ensemble.
    pub isotensional: isotensional::FJC,

    /// The thermodynamic functions of the model in the modified canonical ensemble.
    pub modified_canonical: modified_canonical::FJC
}

/// The implemented functionality of the thermodynamics of the FJC model.
impl FJC
{
    /// Initializes and returns an instance of the thermodynamics of the FJC model.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
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
