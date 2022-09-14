use crate::physics::single_chain::
{
    Isometric,
    IsometricLegendre
};

pub mod test;

pub struct FJC
{
    pub link_length: f64,
    pub number_of_links: u16,
    pub number_of_links_f64: f64,
    pub legendre: FJCLegendre
}

impl Isometric for FJC
{
    fn init(number_of_links: u16, link_length: f64) -> FJC
    {
        FJC
        {
            link_length: link_length,
            number_of_links: number_of_links,
            number_of_links_f64: number_of_links as f64,
            legendre: FJCLegendre::init(number_of_links, link_length)
        }
    }
}

pub struct FJCLegendre
{
    pub link_length: f64,
    pub number_of_links: u16,
    pub number_of_links_f64: f64
}

impl IsometricLegendre for FJCLegendre
{
    fn init(number_of_links: u16, link_length: f64) -> FJCLegendre
    {
        FJCLegendre
        {
            link_length: link_length,
            number_of_links: number_of_links,
            number_of_links_f64: number_of_links as f64
        }
    }
}