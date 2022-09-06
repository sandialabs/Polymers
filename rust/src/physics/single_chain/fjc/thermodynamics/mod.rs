use crate::physics::single_chain::
{
    Math,
    langevin,
    ln_sinhc
};

pub mod test;

pub fn nondimensional_end_to_end_length_per_link<T>(nondimensional_energy: T) -> T
where T:
    Math<T>
{
    langevin(&nondimensional_energy)
}

pub fn nondimensional_gibbs_free_energy_per_link<T>(nondimensional_energy: T) -> T
where T:
    Math<T>,
    T: std::marker::Copy,
    T: std::ops::Mul<T, Output = T>,
    T: std::ops::Sub<T, Output = T>,
{
    langevin(&nondimensional_energy)*nondimensional_energy - ln_sinhc(&nondimensional_energy)
}
