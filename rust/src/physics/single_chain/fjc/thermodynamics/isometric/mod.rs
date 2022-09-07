pub mod test;

pub struct Isometric
{
    pub number_of_links: u16
}

impl Isometric
{
    pub fn init(number_of_links: u16) -> Isometric
    {
        Isometric
        {
            number_of_links: number_of_links
        }
    }
}