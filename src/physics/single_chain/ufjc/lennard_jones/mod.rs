// #[cfg(feature = "python")]
// pub mod py;

mod test;

/// The Lennard-Jones link potential freely-jointed chain (Lennard-Jones-FJC) model thermodynamics.
pub mod thermodynamics;

/// The structure of the Lennard-Jones-FJC model.
pub struct LENNARDJONESFJC
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
    pub thermodynamics: self::thermodynamics::LENNARDJONESFJC
}

/// The implemented functionality of the Lennard-Jones-FJC model.
impl LENNARDJONESFJC
{
    /// Initializes and returns an instance of the Lennard-Jones-FJC model.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        LENNARDJONESFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            thermodynamics: self::thermodynamics::LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness),
        }
    }
}
