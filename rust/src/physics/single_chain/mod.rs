pub mod test;
pub mod fjc;
pub mod swfjc;
use self::{
    fjc::{
        FJC,
        thermodynamics::{
            Isometric as FJCIsometric,
            Isotensional as FJCIsotensional,
            ModifiedCanonical as FJCModifiedCanonical,
            isometric::Legendre as FJCIsometricLegendre,
            isotensional::Legendre as FJCIsotensionalLegendre,
            modified_canonical::
            {
                asymptotic::
                {
                    WeakPotential as ModifiedCanonicalAsymptoticWeakPotential,
                    StrongPotential as ModifiedCanonicalAsymptoticStrongPotential
                }
            }
        }
    },
    swfjc::{
        SWFJC,
        thermodynamics::{
            Isotensional as SWFJCIsotensional,
            isotensional::Legendre as SWFJCIsotensionalLegendre,
        }
    },
};