use polymers::physics::single_chain::fjc::thermodynamics::isometric::
{
    nondimensional_equilibrium_radial_distribution as check_geq,
    monte_carlo::nondimensional_equilibrium_radial_distribution
};

const NUMBER_OF_BINS: usize = 1000;
const NUMBER_OF_LINKS: usize = 8;
const NUMBER_OF_SAMPLES: usize = 10000000;
const TOL: f64 = 5e-2;

#[test]
fn monte_carlo_nondimensional_equilibrium_radial_distribution()
{
    let (gamma, g_eq, _, _, _) = nondimensional_equilibrium_radial_distribution::<NUMBER_OF_BINS, NUMBER_OF_LINKS>(NUMBER_OF_SAMPLES);
    let mut check = 0.0;
    let mut residual = 0.0;
    gamma.iter().zip(g_eq.iter()).for_each(|(gamma_i, g_eq_i)|{
        check = check_geq(&(NUMBER_OF_LINKS as u8), &gamma_i);
        residual = (g_eq_i - &check).abs();
        assert!(residual < TOL || residual/check < TOL || g_eq_i < &TOL);
    });
}

#[test]
fn monte_carlo_isotensional_mechanics_from_first_moments()
{
    let (gamma, g_eq, first_moments, _, _) = nondimensional_equilibrium_radial_distribution::<NUMBER_OF_BINS, NUMBER_OF_LINKS>(NUMBER_OF_SAMPLES);
    // Z = 2*pi*int-dgamma ...
    // <gamma> = sum_i <cos-theta_i>/N_b = ???
    // test small, medium, large force
    todo!()
}
