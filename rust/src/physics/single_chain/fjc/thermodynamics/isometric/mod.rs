use crate::math::
{
    Math,
    langevin
};
use crate::physics::single_chain::Isometric;

pub mod test;

pub struct FJC
{
    pub number_of_links: u16,
    pub link_length: f64
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
    fn force<T>(&self, end_to_end_length: &T) -> T
    where T:
        Math<T>
    {
        langevin(end_to_end_length)
    }
}