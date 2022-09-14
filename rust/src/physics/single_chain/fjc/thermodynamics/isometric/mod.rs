use crate::physics::single_chain::Isometric;

pub mod test;

pub struct FJC
{
    pub link_length: f64,
    pub number_of_links: u16
}

impl Isometric for FJC
{
    fn init(number_of_links: u16, link_length: f64) -> FJC
    {
        FJC
        {
            link_length: link_length,
            number_of_links: number_of_links
        }
    }
}