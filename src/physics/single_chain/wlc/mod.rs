#[cfg(feature = "python")]
pub mod py;

mod test;

/// The worm-like chain (WLC) model thermodynamics.
pub mod thermodynamics;

/// The structure of the WLC model.
pub struct WLC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The persistance length of the chain in units of nm.
    pub persistance_length: f64,

    /// The thermodynamic functions of the model.
    pub thermodynamics: self::thermodynamics::WLC
}

/// The implemented functionality of the WLC model.
impl WLC
{
    /// Initializes and returns an instance of the WLC model.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, persistance_length: f64) -> Self
    {
        WLC
        {
            hinge_mass,
            link_length,
            number_of_links,
            persistance_length,
            thermodynamics: self::thermodynamics::WLC::init(number_of_links, link_length, hinge_mass, persistance_length),
        }
    }
}
