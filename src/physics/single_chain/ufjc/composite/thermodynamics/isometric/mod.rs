#[cfg(feature = "python")]
pub mod py;

mod test;

/// The composite uFJC (CuFJC) model thermodynamics in the isometric ensemble approximated using an asymptotic approach.
pub mod asymptotic;

/// The structure of the thermodynamics of the CuFJC model in the isometric ensemble.
pub struct CUFJC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The number of bonds in each link.
    pub number_of_bonds: u8,

    /// The stiffness of each bond in units of J/(mol⋅nm^2).
    pub bond_stiffness: f64,

    /// The energy of each bond in units of J/mol.
    pub bond_energy: f64,

    /// The scission energy of each bond in units of J/mol.
    pub bond_scission_energy: f64,

    /// The attempt frequency of each bond in units of 1/ns.
    pub bond_attempt_frequency: f64,

    /// The thermodynamic functions of the model in the isometric ensemble approximated using an asymptotic approach.
    pub asymptotic: self::asymptotic::CUFJC

}

/// The implemented functionality of the thermodynamics of the CuFJC model in the isometric ensemble.
impl CUFJC
{
    /// Initializes and returns an instance of the thermodynamics of the CuFJC model in the isometric ensemble.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, number_of_bonds: u8, bond_stiffness: f64, bond_energy: f64, bond_scission_energy: f64, bond_attempt_frequency: f64) -> Self
    {
        CUFJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_bonds,
            bond_stiffness,
            bond_energy,
            bond_scission_energy,
            bond_attempt_frequency,
            asymptotic: self::asymptotic::CUFJC::init(number_of_links, link_length, hinge_mass, number_of_bonds, bond_stiffness, bond_energy, bond_scission_energy, bond_attempt_frequency)
        }
    }
}
