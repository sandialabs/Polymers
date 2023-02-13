#[no_mangle]
pub extern fn physics_single_chain_swfjc_thermodynamics_isotensional_legendre_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, well_width: f64, force: f64, temperature: f64) -> f64
{
    super::helmholtz_free_energy(&number_of_links, &link_length, &hinge_mass, &well_width, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_swfjc_thermodynamics_isotensional_legendre_helmholtz_free_energy_per_link(link_length: f64, hinge_mass: f64, well_width: f64, force: f64, temperature: f64) -> f64
{
    super::helmholtz_free_energy_per_link(&link_length, &hinge_mass, &well_width, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_swfjc_thermodynamics_isotensional_legendre_relative_helmholtz_free_energy(number_of_links: u8, link_length: f64, well_width: f64, force: f64, temperature: f64) -> f64
{
    super::relative_helmholtz_free_energy(&number_of_links, &link_length, &well_width, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_swfjc_thermodynamics_isotensional_legendre_relative_helmholtz_free_energy_per_link(link_length: f64, well_width: f64, force: f64, temperature: f64) -> f64
{
    super::relative_helmholtz_free_energy_per_link(&link_length, &well_width, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_swfjc_thermodynamics_isotensional_legendre_nondimensional_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, well_width: f64, nondimensional_force: f64, temperature: f64) -> f64
{
    super::nondimensional_helmholtz_free_energy(&number_of_links, &link_length, &hinge_mass, &well_width, &nondimensional_force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_swfjc_thermodynamics_isotensional_legendre_nondimensional_helmholtz_free_energy_per_link(link_length: f64, hinge_mass: f64, well_width: f64, nondimensional_force: f64, temperature: f64) -> f64
{
    super::nondimensional_helmholtz_free_energy_per_link(&link_length, &hinge_mass, &well_width, &nondimensional_force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_swfjc_thermodynamics_isotensional_legendre_nondimensional_relative_helmholtz_free_energy(number_of_links: u8, link_length: f64, well_width: f64, nondimensional_force: f64) -> f64
{
    super::nondimensional_relative_helmholtz_free_energy(&number_of_links, &link_length, &well_width, &nondimensional_force)
}
#[no_mangle]
pub extern fn physics_single_chain_swfjc_thermodynamics_isotensional_legendre_nondimensional_relative_helmholtz_free_energy_per_link(link_length: f64, well_width: f64, nondimensional_force: f64) -> f64
{
    super::nondimensional_relative_helmholtz_free_energy_per_link(&link_length, &well_width, &nondimensional_force)
}
