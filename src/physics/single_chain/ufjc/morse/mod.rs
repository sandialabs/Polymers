#[cfg(feature = "python")]
pub mod py;

mod test;

/// The Morse link potential freely-jointed chain (Morse-FJC) model thermodynamics.
pub mod thermodynamics;

/// The structure of the Morse-FJC model.
pub struct MORSEFJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The stiffness of each link in the chain in units of J/(mol⋅nm^2).
    pub link_stiffness: f64,

    /// The energy of each link in the chain in units of J/mol.
    pub link_energy: f64,

    /// The thermodynamic functions of the model.
    pub thermodynamics: self::thermodynamics::MORSEFJC
}

/// The implemented functionality of the Morse-FJC model.
impl MORSEFJC
{
    /// Initializes and returns an instance of the Morse-FJC model.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64, link_energy: f64) -> Self
    {
        MORSEFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            link_energy,
            thermodynamics: self::thermodynamics::MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy),
        }
    }
}
