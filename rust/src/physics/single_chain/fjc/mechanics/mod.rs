use ndarray::Array1;

pub mod test;

pub struct Mechanics
{
    pub reduced_mass: f64,
    pub link_length: f64,
    pub number_of_links: u16,
}

impl Mechanics
{
    pub fn init(number_of_links: u16, link_length: f64, reduced_mass: f64) -> Mechanics
    {
        Mechanics
        {
            reduced_mass,
            link_length,
            number_of_links
        }
    }
    pub fn random_configuration(&self) -> Array1<f64>
    {
        // put random unit vector(s) generator in math, return here, try not to use ndarray
        Array1::<f64>::linspace(1.0, 10.0, (self.number_of_links + 1) as usize)
    }
}

// DOF-specific stuff like potential energy function
// and stuff for Monte Carlo (but keep general stuff in calculation/)
// can keep stuff that's not general but used for certain/all single-chain modesl in single_chain/
