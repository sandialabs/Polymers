#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isometric_asymptotic_legendre_force(number_of_links: u8, link_length: f64, link_stiffness: f64, end_to_end_length: f64, temperature: f64) -> f64
{
    super::force(&number_of_links, &link_length, &link_stiffness, &end_to_end_length, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isometric_asymptotic_legendre_nondimensional_force(nondimensional_link_stiffness: f64, nondimensional_end_to_end_length_per_link: f64) -> f64
{
    super::nondimensional_force(&nondimensional_link_stiffness, &nondimensional_end_to_end_length_per_link)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isometric_asymptotic_legendre_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64, end_to_end_length: f64, temperature: f64) -> f64
{
    super::helmholtz_free_energy(&number_of_links, &link_length, &hinge_mass, &link_stiffness, &end_to_end_length, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isometric_asymptotic_legendre_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64, end_to_end_length: f64, temperature: f64) -> f64
{
    super::helmholtz_free_energy_per_link(&number_of_links, &link_length, &hinge_mass, &link_stiffness, &end_to_end_length, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isometric_asymptotic_legendre_relative_helmholtz_free_energy(number_of_links: u8, link_length: f64, link_stiffness: f64, end_to_end_length: f64, temperature: f64) -> f64
{
    super::relative_helmholtz_free_energy(&number_of_links, &link_length, &link_stiffness, &end_to_end_length, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isometric_asymptotic_legendre_relative_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, link_stiffness: f64, end_to_end_length: f64, temperature: f64) -> f64
{
    super::relative_helmholtz_free_energy_per_link(&number_of_links, &link_length, &link_stiffness, &end_to_end_length, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isometric_asymptotic_legendre_nondimensional_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_link_stiffness: f64, nondimensional_end_to_end_length_per_link: f64, temperature: f64) -> f64
{
    super::nondimensional_helmholtz_free_energy(&number_of_links, &link_length, &hinge_mass, &nondimensional_link_stiffness, &nondimensional_end_to_end_length_per_link, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isometric_asymptotic_legendre_nondimensional_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_link_stiffness: f64, nondimensional_end_to_end_length_per_link: f64, temperature: f64) -> f64
{
    super::nondimensional_helmholtz_free_energy_per_link(&number_of_links, &link_length, &hinge_mass, &nondimensional_link_stiffness, &nondimensional_end_to_end_length_per_link, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isometric_asymptotic_legendre_nondimensional_relative_helmholtz_free_energy(number_of_links: u8, nondimensional_link_stiffness: f64, nondimensional_end_to_end_length_per_link: f64) -> f64
{
    super::nondimensional_relative_helmholtz_free_energy(&number_of_links, &nondimensional_link_stiffness, &nondimensional_end_to_end_length_per_link)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_lennard_jones_thermodynamics_isometric_asymptotic_legendre_nondimensional_relative_helmholtz_free_energy_per_link(nondimensional_link_stiffness: f64, nondimensional_end_to_end_length_per_link: f64) -> f64
{
    super::nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_link_stiffness, &nondimensional_end_to_end_length_per_link)
}