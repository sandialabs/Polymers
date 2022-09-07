use crate::physics::single_chain::
{
    Math,
    langevin,
    ln_sinhc
};

pub mod test;

pub fn nondimensional_end_to_end_length_per_link<T>(nondimensional_force: T) -> T
where T:
    Math<T>
{
    langevin(&nondimensional_force)
}

pub fn nondimensional_relative_gibbs_free_energy_per_link<T>(nondimensional_force: T) -> T
where T:
    Math<T>,
    T: std::marker::Copy,
    T: std::ops::Mul<T, Output = T>,
    T: std::ops::Sub<T, Output = T>,
{
    langevin(&nondimensional_force)*nondimensional_force - ln_sinhc(&nondimensional_force)
}