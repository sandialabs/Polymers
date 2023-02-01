#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isometric_force(number_of_links: u8, link_length: f64, hinge_mass: f64, end_to_end_length: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).force(&end_to_end_length, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isometric_nondimensional_force(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_end_to_end_length_per_link: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_force(&nondimensional_end_to_end_length_per_link)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isometric_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, end_to_end_length: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).helmholtz_free_energy(&end_to_end_length, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isometric_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, end_to_end_length: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).helmholtz_free_energy_per_link(&end_to_end_length, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isometric_relative_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, end_to_end_length: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).relative_helmholtz_free_energy(&end_to_end_length, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isometric_relative_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, end_to_end_length: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isometric_nondimensional_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_end_to_end_length_per_link: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isometric_nondimensional_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_end_to_end_length_per_link: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isometric_nondimensional_relative_helmholtz_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_end_to_end_length_per_link: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isometric_nondimensional_relative_helmholtz_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_end_to_end_length_per_link: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isometric_equilibrium_distribution(number_of_links: u8, link_length: f64, hinge_mass: f64, end_to_end_length: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).equilibrium_distribution(&end_to_end_length)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isometric_nondimensional_equilibrium_distribution(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_end_to_end_length_per_link: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_equilibrium_distribution(&nondimensional_end_to_end_length_per_link)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isometric_equilibrium_radial_distribution(number_of_links: u8, link_length: f64, hinge_mass: f64, end_to_end_length: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).equilibrium_radial_distribution(&end_to_end_length)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isometric_nondimensional_equilibrium_radial_distribution(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_end_to_end_length_per_link: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_equilibrium_radial_distribution(&nondimensional_end_to_end_length_per_link)
}