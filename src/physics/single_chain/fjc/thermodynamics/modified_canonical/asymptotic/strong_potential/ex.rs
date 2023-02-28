#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_modified_canonical_force(number_of_links: u8, link_length: f64, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> f64
{
    super::force(&number_of_links, &link_length, &potential_distance, &potential_stiffness, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_modified_canonical_nondimensional_force(number_of_links: u8, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64) -> f64
{
    super::nondimensional_force(&number_of_links, &nondimensional_potential_distance, &nondimensional_potential_stiffness)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_modified_canonical_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> f64
{
    super::helmholtz_free_energy(&number_of_links, &link_length, &hinge_mass, &potential_distance, &potential_stiffness, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_modified_canonical_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> f64
{
    super::helmholtz_free_energy_per_link(&number_of_links, &link_length, &hinge_mass, &potential_distance, &potential_stiffness, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_modified_canonical_relative_helmholtz_free_energy(number_of_links: u8, link_length: f64, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> f64
{
    super::relative_helmholtz_free_energy(&number_of_links, &link_length, &potential_distance, &potential_stiffness, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_modified_canonical_relative_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, potential_distance: f64, potential_stiffness: f64, temperature: f64) -> f64
{
    super::relative_helmholtz_free_energy_per_link(&number_of_links, &link_length, &potential_distance, &potential_stiffness, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_modified_canonical_nondimensional_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64, temperature: f64) -> f64
{
    super::nondimensional_helmholtz_free_energy(&number_of_links, &link_length, &hinge_mass, &nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_modified_canonical_nondimensional_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64, temperature: f64) -> f64
{
    super::nondimensional_helmholtz_free_energy_per_link(&number_of_links, &link_length, &hinge_mass, &nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_modified_canonical_nondimensional_relative_helmholtz_free_energy(number_of_links: u8, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64) -> f64
{
    super::nondimensional_relative_helmholtz_free_energy(&number_of_links, &nondimensional_potential_distance, &nondimensional_potential_stiffness)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_modified_canonical_nondimensional_relative_helmholtz_free_energy_per_link(number_of_links: u8, nondimensional_potential_distance: f64, nondimensional_potential_stiffness: f64) -> f64
{
    super::nondimensional_relative_helmholtz_free_energy_per_link(&number_of_links, &nondimensional_potential_distance, &nondimensional_potential_stiffness)
}