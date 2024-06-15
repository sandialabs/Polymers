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
    let (gamma, g_eq) = nondimensional_equilibrium_radial_distribution::<1_000, 32>(&(PI/12.0), 2_000_000_000);
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
    let handle_4 = thread::spawn(|| compute(4));
    let handle_5 = thread::spawn(|| compute(5));
    let handle_6 = thread::spawn(|| compute(6));
    let handle_7 = thread::spawn(|| compute(7));
    let handle_8 = thread::spawn(|| compute(8));
    let handle_9 = thread::spawn(|| compute(9));
    let handle_10 = thread::spawn(|| compute(10));
    let handle_11 = thread::spawn(|| compute(11));
    let handle_12 = thread::spawn(|| compute(12));
    let handle_13 = thread::spawn(|| compute(13));
    let handle_14 = thread::spawn(|| compute(14));
    let handle_15 = thread::spawn(|| compute(15));
    let handle_16 = thread::spawn(|| compute(16));
    let handle_17 = thread::spawn(|| compute(17));
    let handle_18 = thread::spawn(|| compute(18));
    let handle_19 = thread::spawn(|| compute(19));
    let handle_20 = thread::spawn(|| compute(20));
    let handle_21 = thread::spawn(|| compute(21));
    let handle_22 = thread::spawn(|| compute(22));
    let handle_23 = thread::spawn(|| compute(23));
    let handle_24 = thread::spawn(|| compute(24));
    compute(25);
    handle_1.join().unwrap();
    handle_2.join().unwrap();
    handle_3.join().unwrap();
    handle_4.join().unwrap();
    handle_5.join().unwrap();
    handle_6.join().unwrap();
    handle_7.join().unwrap();
    handle_8.join().unwrap();
    handle_9.join().unwrap();
    handle_10.join().unwrap();
    handle_11.join().unwrap();
    handle_12.join().unwrap();
    handle_13.join().unwrap();
    handle_14.join().unwrap();
    handle_15.join().unwrap();
    handle_16.join().unwrap();
    handle_17.join().unwrap();
    handle_18.join().unwrap();
    handle_19.join().unwrap();
    handle_20.join().unwrap();
    handle_21.join().unwrap();
    handle_22.join().unwrap();
    handle_23.join().unwrap();
    handle_24.join().unwrap();
}
