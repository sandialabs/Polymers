#[cfg(feature = "python")]
pub mod py;

mod test;

/// The extensible freely-jointed chain (EFJC) model thermodynamics in the isometric ensemble.
pub mod isometric;

/// The extensible freely-jointed chain (EFJC) model thermodynamics in the isotensional ensemble.
pub mod isotensional;

/// The structure of the thermodynamics of the EFJC model.
pub struct EFJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The stiffness of each link in the chain in units of J/(mol⋅nm^2).
    pub link_stiffness: f64,

    /// The thermodynamic functions of the model in the isometric ensemble.
    pub isometric: isometric::EFJC,

    /// The thermodynamic functions of the model in the isotensional ensemble.
    pub isotensional: isotensional::EFJC
}

/// The implemented functionality of the thermodynamics of the EFJC model.
impl EFJC
{
    /// Initializes and returns an instance of the thermodynamics of the EFJC model.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        EFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            isometric: self::isometric::EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness),
            isotensional: self::isotensional::EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness)
        }
    }
}
