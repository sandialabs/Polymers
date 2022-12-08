pub mod test;
pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u8,
}
impl FJC
{
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64) -> FJC
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links
        }
    }
    pub fn random_configuration(&self)
    {
        // put random unit vector(s) generator in math, return here, try not to use ndarray
    }
}
//
// do NOT use rand or other external dependencies
//
// DOF-specific stuff like potential energy function
// and stuff for Monte Carlo (but keep general stuff in calculation/)
// can keep stuff that's not general but used for certain/all single-chain modesl in single_chain/