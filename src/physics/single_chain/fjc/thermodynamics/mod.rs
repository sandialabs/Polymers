#[cfg(feature = "python")]
pub mod py;

mod test;

/// The freely-jointed chain (FJC) model thermodynamics in the isometric ensemble.
pub mod isometric;

/// The freely-jointed chain (FJC) model thermodynamics in the isotensional ensemble.
pub mod isotensional;

/// The freely-jointed chain (FJC) model thermodynamics in the modified canonical ensemble.
pub mod modified_canonical;

/// The structure of the thermodynamics of the FJC model.
pub struct FJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The thermodynamic functions of the model in the isometric ensemble.
    pub isometric: isometric::FJC,

    /// The thermodynamic functions of the model in the isotensional ensemble.
    pub isotensional: isotensional::FJC,

    /// The thermodynamic functions of the model in the modified canonical ensemble.
    pub modified_canonical: modified_canonical::FJC
}

/// The implemented functionality of the thermodynamics of the FJC model.
impl FJC
{
    /// Initializes and returns an instance of the thermodynamics of the FJC model.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> Self
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            isometric: isometric::FJC::init(number_of_links, link_length, hinge_mass),
            isotensional: isotensional::FJC::init(number_of_links, link_length, hinge_mass),
            modified_canonical: modified_canonical::FJC::init(number_of_links, link_length, hinge_mass)
        }
    }
}

fn treloar_sums(number_of_links: &u8, nondimensional_end_to_end_length_per_link: &f64, orders: &Vec<i32>) -> Vec<f64>
{
    let number_of_links_f64 = *number_of_links as f64;
    let n = *number_of_links as u128;
    let p: i32 = (number_of_links - 2).into();
    let m = -*nondimensional_end_to_end_length_per_link*0.5 + 0.5;
    let k = (number_of_links_f64*m).ceil() as u128;
    orders.iter().map(|order| (0..=k-1).collect::<Vec::<u128>>().iter().map(|s| (-1.0_f64).powf(*s as f64)*(((1..=n).product::<u128>()/(1..=*s).product::<u128>()/(1..=n-s).product::<u128>()) as f64)*(m - (*s as f64)/number_of_links_f64).powi(p - order)).sum()).collect()
}