pub mod test;
pub mod fjc;

pub trait Isometric
{
    fn init(number_of_links: u16, link_length: f64, hinge_mass: f64) -> Self;
    fn force(&self, end_to_end_length: &f64, temperature: f64) -> f64;
    fn nondimensional_force(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64;
    fn helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: f64) -> f64;
    fn helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: f64) -> f64;
    fn relative_helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: f64) -> f64;
    fn relative_helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: f64) -> f64;
    fn nondimensional_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: f64) -> f64;
    fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: f64) -> f64;
    fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64;
    fn nondimensional_relative_helmholtz_free_energy_per_link(&self, end_to_end_length: &f64) -> f64;
    fn equilibrium_distribution(&self, end_to_end_length: &f64) -> f64;
    fn nondimensional_equilibrium_distribution(&self, end_to_end_length: &f64) -> f64;
    fn equilibrium_radial_distribution(&self, end_to_end_length: &f64) -> f64;
    fn nondimensional_equilibrium_radial_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64;
}

pub trait Isotensional
{
    fn init(number_of_links: u16, link_length: f64, hinge_mass: f64) -> Self;
    fn end_to_end_length(&self, force: &f64, temperature: f64) -> f64;
    fn end_to_end_length_per_link(&self, force: &f64, temperature: f64) -> f64;
    fn nondimensional_end_to_end_length(&self, nondimensional_force: &f64) -> f64;
    fn nondimensional_end_to_end_length_per_link(&self, nondimensional_force: &f64) -> f64;
    fn gibbs_free_energy(&self, force: &f64, temperature: f64) -> f64;
    fn gibbs_free_energy_per_link(&self, force: &f64, temperature: f64) -> f64;
    fn relative_gibbs_free_energy(&self, force: &f64, temperature: f64) -> f64;
    fn relative_gibbs_free_energy_per_link(&self, force: &f64, temperature: f64) -> f64;
    fn nondimensional_gibbs_free_energy(&self, nondimensional_force: &f64, temperature: f64) -> f64;
    fn nondimensional_gibbs_free_energy_per_link(&self, nondimensional_force: &f64, temperature: f64) -> f64;
    fn nondimensional_relative_gibbs_free_energy(&self, nondimensional_force: &f64) -> f64;
    fn nondimensional_relative_gibbs_free_energy_per_link(&self, nondimensional_force: &f64) -> f64;
}

pub trait IsometricLegendre
{
    fn init(number_of_links: u16, link_length: f64, hinge_mass: f64) -> Self;
    fn force(&self, end_to_end_length: &f64, temperature: f64) -> f64;
    fn nondimensional_force(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64;
    fn helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: f64) -> f64;
    fn helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: f64) -> f64;
    fn relative_helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: f64) -> f64;
    fn relative_helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: f64) -> f64;
    fn nondimensional_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: f64) -> f64;
    fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: f64) -> f64;
    fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64;
    fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64;
}

pub trait IsotensionalLegendre
{
    fn init(number_of_links: u16, link_length: f64, hinge_mass: f64) -> Self;
    fn helmholtz_free_energy(&self, force: &f64, temperature: f64) -> f64;
    fn helmholtz_free_energy_per_link(&self, force: &f64, temperature: f64) -> f64;
    fn relative_helmholtz_free_energy(&self, force: &f64, temperature: f64) -> f64;
    fn relative_helmholtz_free_energy_per_link(&self, force: &f64, temperature: f64) -> f64;
    fn nondimensional_helmholtz_free_energy(&self, nondimensional_force: &f64, temperature: f64) -> f64;
    fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64, temperature: f64) -> f64;
    fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_force: &f64) -> f64;
    fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_force: &f64) -> f64;
}
