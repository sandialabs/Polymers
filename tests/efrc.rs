use polymers::physics::single_chain::
{
    frc::thermodynamics::isometric::monte_carlo::nondimensional_equilibrium_radial_distribution as frc_nondimensional_equilibrium_radial_distribution,
    efrc::thermodynamics::isometric::monte_carlo::nondimensional_equilibrium_radial_distribution
};

use std::f64::consts::PI;

const KAPPA_LARGE: f64 = 1e5;
const NUMBER_OF_BINS: usize = 1_000;
const NUMBER_OF_LINKS: usize = 8;
const NUMBER_OF_SAMPLES: usize = 10_000_000;
const THETA: f64 = PI/4.0;
const TOL: f64 = 5e-2;

#[test]
fn monte_carlo_nondimensional_equilibrium_radial_distribution_frc_limit()
{
    let gamma_max = 1.01 * (2.0 - 2.0*(PI - THETA).cos()).sqrt()/2.0;
    let (gamma, g_eq) = nondimensional_equilibrium_radial_distribution::<NUMBER_OF_BINS, NUMBER_OF_LINKS>(&gamma_max, &KAPPA_LARGE, &THETA, NUMBER_OF_SAMPLES);
    let (gamma_frc, g_eq_frc) = frc_nondimensional_equilibrium_radial_distribution::<NUMBER_OF_BINS, NUMBER_OF_LINKS>(&THETA, NUMBER_OF_SAMPLES);
    let mut residual = 0.0;
    gamma_frc.iter().zip(g_eq_frc.iter()).for_each(|(gamma_frc_i, g_eq_frc_i)|{
        for j in 1..NUMBER_OF_BINS {
            if gamma_frc_i < &gamma[j] {
                residual = (g_eq[j - 1] + (g_eq[j] - g_eq[j - 1])/(gamma[j] - gamma[j - 1])*(gamma_frc_i - gamma[j - 1]) - g_eq_frc_i).abs();
                assert!(residual < TOL || residual/g_eq_frc_i < TOL);
                break
            }
        }
    });
}
