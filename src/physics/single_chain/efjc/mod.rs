#[cfg(feature = "python")]
pub mod py;

mod test;

/// The extensible freely-jointed chain (EFJC) model thermodynamics.
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

    /// The stiffness of each link in the chain in units of J/(mol⋅nm^2).
    pub link_stiffness: f64,

    /// The thermodynamic functions of the model.
    pub thermodynamics: self::thermodynamics::EFJC
}

/// The implemented functionality of the EFJC model.
impl EFJC
{
    /// Initializes and returns an instance of the EFJC model.
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
