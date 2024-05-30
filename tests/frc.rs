use polymers::physics::single_chain::
{
    frc::thermodynamics::isometric::monte_carlo::nondimensional_equilibrium_radial_distribution,
    wlc::thermodynamics::isometric::WLC
};

const KAPPA: f64 = 5.0/11.0;
const NUMBER_OF_BINS: usize = 1000;
const NUMBER_OF_LINKS: usize = 256;
const NUMBER_OF_SAMPLES: usize = 10000000;
const TOL: f64 = 8e-2;

#[test]
fn temporary()
{
    let (gamma, g_eq) = nondimensional_equilibrium_radial_distribution::<1_000, 8>(&(std::f64::consts::PI/3.0), 1_000_000_000);
    gamma.iter().zip(g_eq.iter()).for_each(|(gamma, g_eq)|
        println!("{}\t{}", gamma, g_eq)
    );
}

#[test]
fn todo()
{
    todo!("WLC analytic implementation seems to be wrong w.r.t. paper it came from and the Monte Carlo results.")
}

#[test]
fn monte_carlo_nondimensional_equilibrium_radial_distribution_wlc_limit()
{
    let theta = (2.0/(NUMBER_OF_LINKS as f64)/KAPPA).sqrt();
    let (gamma, g_eq) = nondimensional_equilibrium_radial_distribution::<NUMBER_OF_BINS, NUMBER_OF_LINKS>(&theta, NUMBER_OF_SAMPLES);
    let mut check = 0.0;
    let mut residual = 0.0;
    let wlc = WLC::init(1, 1.0, 1.0, KAPPA);
    gamma.iter().zip(g_eq.iter()).for_each(|(gamma_i, g_eq_i)|{
        check = wlc.nondimensional_equilibrium_radial_distribution(&gamma_i);
        residual = (g_eq_i - &check).abs();
        assert!(residual < TOL || residual/check < TOL || g_eq_i < &TOL);
    });
}
