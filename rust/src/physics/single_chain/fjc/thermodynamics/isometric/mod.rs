pub mod test;

pub struct Isometric
{
    pub number_of_links: u16,
    pub link_length: f64
}

impl Isometric
{
    pub fn init(number_of_links: u16, link_length: f64) -> Isometric
    {
        Isometric
        {
            number_of_links: number_of_links,
            link_length: link_length,
        }
    }
}
