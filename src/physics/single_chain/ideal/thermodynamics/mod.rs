#[cfg(feature = "python")]
pub mod py;

mod test;

/// The ideal chain model thermodynamics in the isometric ensemble.
pub mod isometric;

/// The ideal chain model thermodynamics in the isotensional ensemble.
pub mod isotensional;

/// The structure of the thermodynamics of the ideal chain model.
pub struct Ideal
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The thermodynamic functions of the model in the isometric ensemble.
    pub isometric: isometric::Ideal,

    /// The thermodynamic functions of the model in the isotensional ensemble.
    pub isotensional: isotensional::Ideal
}

/// The implemented functionality of the thermodynamics of the ideal chain model.
impl Ideal
{
    /// Initializes and returns an instance of the thermodynamics of the ideal chain model.
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
