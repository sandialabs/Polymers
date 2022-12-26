pub mod test;
pub mod ideal;
pub mod fjc;
pub mod efjc;
pub mod swfjc;
use super::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT,
    ZERO,
    POINTS
};
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
            isotensional::{
                Legendre as EFJCIsotensionalLegendre,
                Asymptotic as EFJCIsotensionalAsymptotic,
                asymptotic::{
                    Alternative as EFJCIsotensionalAsymptoticAlternative,
                    Reduced as EFJCIsotensionalAsymptoticReduced,
                    Legendre as EFJCIsotensionalAsymptoticLegendre,
                    alternative::Legendre as EFJCIsotensionalAsymptoticAlternativeLegendre,
                    reduced::Legendre as EFJCIsotensionalAsymptoticReducedLegendre
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
