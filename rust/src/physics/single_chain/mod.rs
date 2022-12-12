pub mod test;
pub mod fjc;
pub mod swfjc;
use self::{
    fjc::{
        thermodynamics::{
            Isometric as FJCIsometric,
            Isotensional as FJCIsotensional,
            ModifiedCanonical as FJCModifiedCanonical,
            isometric::Legendre as FJCIsometricLegendre,
            isotensional::Legendre as FJCIsotensionalLegendre,
            modified_canonical::asymptotic::{
                WeakPotential as FJCModifiedCanonicalAsymptoticWeakPotential,
                StrongPotential as FJCModifiedCanonicalAsymptoticStrongPotential
            }
        }
    },
    swfjc::{
        thermodynamics::{
            Isometric as SWFJCIsometric,
            Isotensional as SWFJCIsotensional,
            isometric::Legendre as SWFJCIsometricLegendre,
            isotensional::Legendre as SWFJCIsotensionalLegendre,
        }
    },
};