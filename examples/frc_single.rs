use polymers::physics::single_chain::frc::thermodynamics::isometric::monte_carlo::nondimensional_equilibrium_radial_distribution;

use std::
{
    env,
    f64::consts::PI,
    fs::File,
    io::Write,
    thread
};

fn compute(index: usize)
{
    let file_path = &env::args().collect::<Vec<String>>()[1];
    let (gamma, g_eq) = nondimensional_equilibrium_radial_distribution::<1_000, 8>(&(PI/4.0), 20_000_000);
    let file_name = "data/frc/".to_owned() + &file_path + &format!("/frc_{}.csv", index);
    let mut file = File::create(file_name).unwrap();
    gamma.iter().zip(g_eq.iter()).for_each(|(gamma, g_eq)|{
        writeln!(&mut file, "{}\t{}", gamma, g_eq).unwrap();
    });
}

fn main() {
    compute(0);
}
