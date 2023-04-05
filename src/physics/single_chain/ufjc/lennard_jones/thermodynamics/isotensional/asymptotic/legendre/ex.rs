#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isotensional_asymptotic_legendre_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64, force: f64, temperature: f64) -> f64
{
    super::helmholtz_free_energy(&number_of_links, &link_length, &hinge_mass, &link_stiffness, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isotensional_asymptotic_legendre_helmholtz_free_energy_per_link(link_length: f64, hinge_mass: f64, link_stiffness: f64, force: f64, temperature: f64) -> f64
{
    super::helmholtz_free_energy_per_link(&link_length, &hinge_mass, &link_stiffness, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isotensional_asymptotic_legendre_relative_helmholtz_free_energy(number_of_links: u8, link_length: f64, link_stiffness: f64, force: f64, temperature: f64) -> f64
{
    super::relative_helmholtz_free_energy(&number_of_links, &link_length, &link_stiffness, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isotensional_asymptotic_legendre_relative_helmholtz_free_energy_per_link(link_length: f64, link_stiffness: f64, force: f64, temperature: f64) -> f64
{
    super::relative_helmholtz_free_energy_per_link(&link_length, &link_stiffness, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isotensional_asymptotic_legendre_nondimensional_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_link_stiffness: f64, nondimensional_force: f64, temperature: f64) -> f64
{
    super::nondimensional_helmholtz_free_energy(&number_of_links, &link_length, &hinge_mass, &nondimensional_link_stiffness, &nondimensional_force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isotensional_asymptotic_legendre_nondimensional_helmholtz_free_energy_per_link(link_length: f64, hinge_mass: f64, nondimensional_link_stiffness: f64, nondimensional_force: f64, temperature: f64) -> f64
{
    super::nondimensional_helmholtz_free_energy_per_link(&link_length, &hinge_mass, &nondimensional_link_stiffness, &nondimensional_force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isotensional_asymptotic_legendre_nondimensional_relative_helmholtz_free_energy(number_of_links: u8, nondimensional_link_stiffness: f64, nondimensional_force: f64) -> f64
{
    super::nondimensional_relative_helmholtz_free_energy(&number_of_links, &nondimensional_link_stiffness, &nondimensional_force)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isotensional_asymptotic_legendre_nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_link_stiffness: f64, nondimensional_force: f64) -> f64
{
    super::nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_link_stiffness, &nondimensional_force)
}