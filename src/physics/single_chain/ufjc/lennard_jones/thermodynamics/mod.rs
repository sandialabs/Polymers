#[cfg(feature = "python")]
pub mod py;

mod test;

/// The Lennard-Jones link potential freely-jointed chain (Lennard-Jones-FJC) model thermodynamics in the isotensional ensemble.
pub mod isotensional;

/// The structure of the Lennard-Jones-FJC model thermodynamics.
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

    /// The thermodynamic functions of the model in the isotensional ensemble.
    pub isotensional: self::isotensional::LENNARDJONESFJC
}

/// The implemented functionality of the Lennard-Jones-FJC model thermodynamics.
impl LENNARDJONESFJC
{
    /// Initializes and returns an instance of the Lennard-Jones-FJC model thermodynamics.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        LENNARDJONESFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            isotensional: self::isotensional::LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness),
        }
    }
}

fn nondimensional_bond_stretch(nondimensional_link_stiffness: &f64, nondimensional_force: &f64) -> f64
{
    let p = -1.0/13.0;
    let s = 6.0*nondimensional_force/nondimensional_link_stiffness;
    let mut lambda: f64 = 1.1;
    let mut residual_rel = 1.0;
    while residual_rel > 1e-5
    {
        lambda = (lambda.powi(-7) - s).powf(p);
        residual_rel = (1.0 - (lambda.powi(-7) - lambda.powi(-13))/s).abs();
    }
    lambda
}