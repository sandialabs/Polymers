use ndarray::Array1;

pub mod test;

pub struct Mechanics
{
    pub number_of_links: u16
}

impl Mechanics
{
    pub fn init(number_of_links: u16) -> Mechanics
    {
        Mechanics
        {
            number_of_links: number_of_links
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