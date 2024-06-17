use polymers::physics::single_chain::
{
    fjc::thermodynamics::isometric::monte_carlo::nondimensional_equilibrium_radial_distribution as frc_nondimensional_equilibrium_radial_distribution,
    efjc::thermodynamics::isometric::monte_carlo::nondimensional_equilibrium_radial_distribution
};

const GAMMA_MAX: f64 = 1.01;
const KAPPA_LARGE: f64 = 1e5;
const NUMBER_OF_BINS: usize = 1_000;
const NUMBER_OF_LINKS: usize = 8;
const NUMBER_OF_SAMPLES: usize = 10_000_000;
const TOL: f64 = 5e-2;

#[test]
fn monte_carlo_nondimensional_equilibrium_radial_distribution_fjc_limit()
{
    let (gamma, g_eq) = nondimensional_equilibrium_radial_distribution::<NUMBER_OF_BINS, NUMBER_OF_LINKS>(&GAMMA_MAX, &KAPPA_LARGE, NUMBER_OF_SAMPLES);
    let (gamma_frc, g_eq_frc) = frc_nondimensional_equilibrium_radial_distribution::<NUMBER_OF_BINS, NUMBER_OF_LINKS>(NUMBER_OF_SAMPLES);
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
