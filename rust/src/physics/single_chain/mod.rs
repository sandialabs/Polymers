use crate::physics::Physics;

pub mod test;
pub mod fjc;

pub struct SingleChain
{
    pub physics: Physics
}

impl SingleChain
{
    pub fn init() -> SingleChain
    {
        SingleChain
        {
            physics: Physics::init()
        }
    }
}