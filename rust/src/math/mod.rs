pub mod test;

pub fn sequential_uniform_random_integer(mut n: u64) -> u64
{
    n = n^(n << 12);
    n = n^(n >> 25);
    n = n^(n << 27);
    n
}

pub fn sequential_uniform_random(mut x: f64) -> f64
{
    x = (sequential_uniform_random_integer(((u64::MAX as f64)*x).ceil() as u64) as f64)/(u64::MAX as f64);
    x
}

pub fn invert<F>(value: f64, function: F, mut guess: f64) -> f64
where F: Fn(f64) -> f64
{
    let mut residual: f64 = 1.0;
    while residual.abs()/value > 1e-4
    {
        residual = function(guess) - value;
        guess = guess - residual/(function(guess + 5e-4) - function(guess - 5e-4))*1e-3;
    }
    guess
}

pub fn integrate<F>(function: F, lower_lim: &f64, upper_lim: &f64, num_points: &u128) -> f64
where F: Fn(f64) -> f64
{
    let dx = (upper_lim - lower_lim)/(*num_points as f64);
    (0..=num_points-1).collect::<Vec::<u128>>().iter().map(|index| function(lower_lim + (0.5 + *index as f64)*dx)).sum::<f64>()*dx
}

pub fn approximate_inverse_langevin(x: &f64) -> f64
{
    (2.14234*x.powf(3.0) - 4.22785*x.powf(2.0) + 3.0*x)/(1.0 - x)/(0.71716*x.powf(3.0) - 0.41103*x.powf(2.0) - 0.39165*x + 1.0)
}

pub fn inverse_langevin(y: &f64, tol: f64) -> f64
{
    let mut x = approximate_inverse_langevin(y);
    let mut residual_rel: f64 = 1.0;
    while residual_rel.abs() > tol
    {
        x = 1.0/(1.0/x.tanh() - y);
        residual_rel = 1.0/x.tanh() - 1.0/x - 1.0;
    }
    x
}
