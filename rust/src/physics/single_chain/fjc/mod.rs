use crate::physics::single_chain::fjc::mechanics::Mechanics;
use crate::physics::single_chain::fjc::thermodynamics::Thermodynamics;

pub mod test;
pub mod mechanics;
pub mod thermodynamics;

pub struct FJC
{
    pub number_of_links: u16,
    pub mechanics: Mechanics,
    pub thermodynamics: Thermodynamics
}

impl FJC
{
    pub fn init(number_of_links: u16) -> FJC
    {
        FJC
        {
            number_of_links: number_of_links,
            mechanics: Mechanics::init(number_of_links),
            thermodynamics: Thermodynamics::init(number_of_links),
        }
    }
}