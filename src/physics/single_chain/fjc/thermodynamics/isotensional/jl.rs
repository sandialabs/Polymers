#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_end_to_end_length(number_of_links: u8, link_length: f64, hinge_mass: f64, force: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).end_to_end_length(&force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_end_to_end_length_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, force: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).end_to_end_length_per_link(&force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_nondimensional_end_to_end_length(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_force: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_end_to_end_length(&nondimensional_force)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_nondimensional_end_to_end_length_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_force: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_end_to_end_length_per_link(&nondimensional_force)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_gibbs_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, force: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).gibbs_free_energy(&force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_gibbs_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, force: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).gibbs_free_energy_per_link(&force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_relative_gibbs_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, force: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).relative_gibbs_free_energy(&force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_relative_gibbs_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, force: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).relative_gibbs_free_energy_per_link(&force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_nondimensional_gibbs_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_force: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_gibbs_free_energy(&nondimensional_force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_nondimensional_gibbs_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_force: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_nondimensional_relative_gibbs_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_force: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_relative_gibbs_free_energy(&nondimensional_force)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_nondimensional_relative_gibbs_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_force: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force)
}