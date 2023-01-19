#[cfg(feature = "python")]
pub mod py;

mod test;

/// The ideal chain model thermodynamics.
pub mod thermodynamics;

/// The structure of the ideal chain model.
pub struct Ideal
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The thermodynamic functions of the model.
    pub thermodynamics: thermodynamics::Ideal
}

/// The implemented functionality of the ideal chain model.
impl Ideal
{
    /// Initializes and returns an instance of the ideal chain model.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        Ideal
        {
            hinge_mass,
            link_length,
            number_of_links,
            thermodynamics: thermodynamics::Ideal::init(number_of_links, link_length, hinge_mass),
        }
    }
}
