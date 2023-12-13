use polymers::physics::single_chain::efjc::thermodynamics::isometric::monte_carlo::nondimensional_equilibrium_radial_distribution;

const GAMMA_MAX: f64 = 1.5;
const KAPPA: f64 = 50.0;
const NUMBER_OF_BINS: usize = 1000;
const NUMBER_OF_LINKS: usize = 3;
const NUMBER_OF_SAMPLES: usize = 10000000;

#[test]
fn monte_carlo_nondimensional_equilibrium_radial_distribution()
{
    let (gamma, g_eq) = nondimensional_equilibrium_radial_distribution::<NUMBER_OF_BINS, NUMBER_OF_LINKS>(&GAMMA_MAX, &KAPPA, NUMBER_OF_SAMPLES);
    gamma.iter().zip(g_eq.iter()).for_each(|output|
        println!("{:?}", output)
    );
}
