pub mod test;
pub mod ideal;
pub mod fjc;
pub mod efjc;
pub mod swfjc;
use self::{
    ideal::{
        Ideal,
        thermodynamics::{
            Isometric as IdealIsometric,
            Isotensional as IdealIsotensional
        }
    },
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
    efjc::{
        EFJC,
        thermodynamics::{
            Isotensional as EFJCIsotensional,
            isotensional::Legendre as EFJCIsotensionalLegendre,
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
