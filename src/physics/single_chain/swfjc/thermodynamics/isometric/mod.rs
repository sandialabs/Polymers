#[cfg(feature = "python")]
pub mod py;

mod test;

/// The square-well freely-jointed chain (SWFJC) model thermodynamics in the isometric ensemble approximated using a Legendre transformation.
pub mod legendre;

/// The structure of the thermodynamics of the SWFJC model in the isometric ensemble.
pub struct SWFJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The width of the well in units of nm.
    pub well_width: f64,

    /// The thermodynamic functions of the model in the isometric ensemble approximated using a Legendre transformation.
    pub legendre: legendre::SWFJC
}

/// The implemented functionality of the SWFJC model in the isometric ensemble.
impl SWFJC
{
    /// Initializes and returns an instance of the SWFJC model in the isometric ensemble.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, well_width: f64) -> Self
    {
        SWFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            well_width,
            legendre: self::legendre::SWFJC::init(number_of_links, link_length, hinge_mass, well_width)
        }
    }
}
