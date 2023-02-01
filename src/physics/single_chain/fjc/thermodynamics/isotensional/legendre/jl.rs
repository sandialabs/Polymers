#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_legendre_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, force: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).helmholtz_free_energy(&force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_legendre_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, force: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).helmholtz_free_energy_per_link(&force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_legendre_relative_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, force: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).relative_helmholtz_free_energy(&force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_legendre_relative_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, force: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).relative_helmholtz_free_energy_per_link(&force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_legendre_nondimensional_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_force: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_helmholtz_free_energy(&nondimensional_force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_legendre_nondimensional_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_force: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_helmholtz_free_energy_per_link(&nondimensional_force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_legendre_nondimensional_relative_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_force: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_relative_helmholtz_free_energy(&nondimensional_force)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_legendre_nondimensional_relative_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_force: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force)
}