use crate::constitutive::Constitutive;

pub mod test;
pub mod arruda_boyce;
pub mod buche_silberstein;
pub mod gent;
pub mod neo_hookean;
pub mod yeoh;

// for more see https://doi.org/10.1016/j.jmbbm.2022.105522

pub struct Hyperelastic {
    pub constitutive: Constitutive
}

impl Hyperelastic {
    pub fn init() -> Hyperelastic {
        Hyperelastic {
            constitutive: Constitutive::init()
        }
    }
}
