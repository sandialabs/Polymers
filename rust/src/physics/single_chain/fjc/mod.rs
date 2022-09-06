use crate::physics::single_chain::SingleChain;

pub mod test;
pub mod thermodynamics;

pub struct FJC
{
    pub single_chain: SingleChain
}

impl FJC
{
    pub fn init() -> FJC
    {
        FJC
        {
            single_chain: SingleChain::init()
        }
    }
}