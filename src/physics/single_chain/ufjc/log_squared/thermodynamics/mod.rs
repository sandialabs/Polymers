#[cfg(feature = "python")]
pub mod py;

mod test;

/// The log-squared link potential freely-jointed chain (log-squared-FJC) model thermodynamics in the isotensional ensemble.
pub mod isotensional;

/// The structure of the log-squared-FJC model thermodynamics.
pub struct LOGSQUAREDFJC
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
    pub isotensional: self::isotensional::LOGSQUAREDFJC
}

/// The implemented functionality of the log-squared-FJC model thermodynamics.
impl LOGSQUAREDFJC
{
    /// Initializes and returns an instance of the log-squared-FJC model thermodynamics.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64) -> Self
    {
        LOGSQUAREDFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            link_stiffness,
            isotensional: self::isotensional::LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness),
        }
    }
}

fn lambert_w(x: &f64) -> f64
{
    let amount_of_iterations = ((x.log10()/3.0).ceil() as u8).max(4_u8);
    let mut w = 0.75*(x + 1.0).ln();
    for _ in 0..amount_of_iterations
    {
        w -= (w*w.exp() - x)/(w.exp()*(w + 1.0) - (w + 2.0)*(w*w.exp() - x)/(2.0 * w + 2.0));
    }
    w
}