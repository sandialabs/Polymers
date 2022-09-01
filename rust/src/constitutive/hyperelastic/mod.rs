use crate::constitutive::Constitutive;

pub mod tests;
pub mod arruda_boyce;
pub mod buche_silberstein;
pub mod gent;
pub mod neo_hookean;
pub mod yeoh;

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
