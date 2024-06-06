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
    let (gamma, g_eq) = nondimensional_equilibrium_radial_distribution::<1_000, 8>(&(PI/4.0), 10_000_000);
    let file_name = file_path.to_owned() + &format!("/frc_{}.csv", index);
    let mut file = File::create(file_name).unwrap();
    gamma.iter().zip(g_eq.iter()).for_each(|(gamma, g_eq)|{
        writeln!(&mut file, "{}\t{}", gamma, g_eq).unwrap();
    });
}

fn main() {
    let handle_1 = thread::spawn(|| compute(1));
    let handle_2 = thread::spawn(|| compute(2));
    let handle_3 = thread::spawn(|| compute(3));
    compute(0);
    handle_1.join().unwrap();
    handle_2.join().unwrap();
    handle_3.join().unwrap();
}