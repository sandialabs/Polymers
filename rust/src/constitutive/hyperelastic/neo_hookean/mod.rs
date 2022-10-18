use crate::constitutive::hyperelastic::Hyperelastic;

pub mod test;

pub struct NeoHookean {
    pub hyperelastic: Hyperelastic
}

impl NeoHookean {
    pub fn init() -> NeoHookean {
        NeoHookean {
            hyperelastic: Hyperelastic::init()
        }
    }
}
