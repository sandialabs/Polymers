pub fn integrate<F>(function: F, lower_lim: &f64, upper_lim: &f64, num_points: &u128) -> f64
where F: Fn(f64) -> f64
{
    let dx = (upper_lim - lower_lim)/(*num_points as f64);
    (0..=num_points-1).collect::<Vec::<u128>>().iter().map(|index| function(lower_lim + (0.5 + *index as f64)*dx)).sum::<f64>()*dx
}
