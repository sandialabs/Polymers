use polymers::physics::single_chain::efrc::thermodynamics::isometric::monte_carlo::nondimensional_equilibrium_radial_distribution;

use std::
{
    f64::consts::PI,
    fs::File,
    io::Write
};

fn main()
{
    let (gamma, g_eq) = nondimensional_equilibrium_radial_distribution::<1_000, 8>(&1.2, &50.0, &(PI/4.0), 1_000_000);
    let mut file = File::create("examples/efrc.csv").unwrap();
    gamma.iter().zip(g_eq.iter()).for_each(|(gamma, g_eq)|{
        writeln!(&mut file, "{}\t{}", gamma, g_eq).unwrap();
    });
}