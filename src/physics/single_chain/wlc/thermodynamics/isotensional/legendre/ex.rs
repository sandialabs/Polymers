#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isotensional_legendre_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, persistance_length: f64, force: f64, temperature: f64) -> f64
{
    super::helmholtz_free_energy(&number_of_links, &link_length, &hinge_mass, &persistance_length, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isotensional_legendre_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, persistance_length: f64, force: f64, temperature: f64) -> f64
{
    super::helmholtz_free_energy_per_link(&number_of_links, &link_length, &hinge_mass, &persistance_length, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isotensional_legendre_relative_helmholtz_free_energy(number_of_links: u8, link_length: f64, persistance_length: f64, force: f64, temperature: f64) -> f64
{
    super::relative_helmholtz_free_energy(&number_of_links, &link_length, &persistance_length, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isotensional_legendre_relative_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, persistance_length: f64, force: f64, temperature: f64) -> f64
{
    super::relative_helmholtz_free_energy_per_link(&number_of_links, &link_length, &persistance_length, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isotensional_legendre_nondimensional_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_persistance_length: f64, nondimensional_force: f64, temperature: f64) -> f64
{
    super::nondimensional_helmholtz_free_energy(&number_of_links, &link_length, &hinge_mass, &nondimensional_persistance_length, &nondimensional_force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isotensional_legendre_nondimensional_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_persistance_length: f64, nondimensional_force: f64, temperature: f64) -> f64
{
    super::nondimensional_helmholtz_free_energy_per_link(&number_of_links, &link_length, &hinge_mass, &nondimensional_persistance_length, &nondimensional_force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isotensional_legendre_nondimensional_relative_helmholtz_free_energy(number_of_links: u8, nondimensional_persistance_length: f64, nondimensional_force: f64) -> f64
{
    super::nondimensional_relative_helmholtz_free_energy(&number_of_links, &nondimensional_persistance_length, &nondimensional_force)
}
#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isotensional_legendre_nondimensional_relative_helmholtz_free_energy_per_link(number_of_links: u8, nondimensional_persistance_length: f64, nondimensional_force: f64) -> f64
{
    super::nondimensional_relative_helmholtz_free_energy_per_link(&number_of_links, &nondimensional_persistance_length, &nondimensional_force)
}