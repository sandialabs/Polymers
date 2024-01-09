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
    for the corrections (not relevant here)
    dont you get terms from derivatives of the sin(theta)s in the integrand jacobian?
    or are those gone somehow?

    let (gamma, g_eq, first_moments, _, _) = nondimensional_equilibrium_radial_distribution::<NUMBER_OF_BINS, NUMBER_OF_LINKS>(NUMBER_OF_SAMPLES);
    let dgamma = gamma[1] - gamma[0];
    let eta = 1.0_f64;
    let number_of_links = NUMBER_OF_LINKS as f64;
    let z = gamma.iter().zip(g_eq.iter()).map(|(gamma_i, g_eq_i)| (number_of_links*eta*gamma_i).sinh()/(number_of_links*eta*gamma_i)*g_eq_i).sum::<f64>()*dgamma;
    first_moments.iter().for_each(|first_moment_j| first_moment_j.iter().for_each(|first_moment_j_i| println!("{}", first_moment_j_i)));
    let n = first_moments.iter().map(|first_moment_j| gamma.iter().zip(g_eq.iter().zip(first_moment_j.iter())).map(|(gamma_i, (g_eq_i, first_moment_j_i))| first_moment_j_i*(number_of_links*eta*gamma_i).sinh()/(number_of_links*eta*gamma_i)*g_eq_i).sum::<f64>()).sum::<f64>()*dgamma;
    println!("{:?}", (z, n, n/z));
    // test small, medium, large force
    todo!()
}
