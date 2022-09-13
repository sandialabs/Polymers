pub mod test;

pub struct Isometric
{
    pub number_of_links: u16,
    pub link_length_in_meters: f64
}

impl Isometric
{
    pub fn init(number_of_links: u16) -> Isometric
    {
        Isometric
        {
            number_of_links: number_of_links,
            link_length_in_meters: 1.0,
        }
    }
}
