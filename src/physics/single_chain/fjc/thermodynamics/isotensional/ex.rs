#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_end_to_end_length(number_of_links: u8, link_length: f64, force: f64, temperature: f64) -> f64
{
    super::end_to_end_length(&number_of_links, &link_length, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_end_to_end_length_per_link(link_length: f64, force: f64, temperature: f64) -> f64
{
    super::end_to_end_length_per_link(&link_length, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_nondimensional_end_to_end_length(number_of_links: u8, nondimensional_force: f64) -> f64
{
    super::nondimensional_end_to_end_length(&number_of_links, &nondimensional_force)
}
#[no_mangle]
pub extern fn physics_single_chain_fjc_thermodynamics_isotensional_nondimensional_end_to_end_length_per_link(nondimensional_force: f64) -> f64
{
    super::nondimensional_end_to_end_length_per_link(&nondimensional_force)
}
