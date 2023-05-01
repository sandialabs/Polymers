#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isometric_gibbs_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, persistance_length: f64, end_to_end_length: f64, temperature: f64) -> f64
{
    super::gibbs_free_energy(&number_of_links, &link_length, &hinge_mass, &persistance_length, &end_to_end_length, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isometric_gibbs_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, persistance_length: f64, end_to_end_length: f64, temperature: f64) -> f64
{
    super::gibbs_free_energy_per_link(&number_of_links, &link_length, &hinge_mass, &persistance_length, &end_to_end_length, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isometric_relative_gibbs_free_energy(number_of_links: u8, link_length: f64, persistance_length: f64, end_to_end_length: f64, temperature: f64) -> f64
{
    super::relative_gibbs_free_energy(&number_of_links, &link_length, &persistance_length, &end_to_end_length, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isometric_relative_gibbs_free_energy_per_link(number_of_links: u8, link_length: f64, persistance_length: f64, end_to_end_length: f64, temperature: f64) -> f64
{
    super::relative_gibbs_free_energy_per_link(&number_of_links, &link_length, &persistance_length, &end_to_end_length, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isometric_nondimensional_gibbs_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_persistance_length: f64, nondimensional_end_to_end_length_per_link: f64, temperature: f64) -> f64
{
    super::nondimensional_gibbs_free_energy(&number_of_links, &link_length, &hinge_mass, &nondimensional_persistance_length, &nondimensional_end_to_end_length_per_link, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isometric_nondimensional_gibbs_free_energy_per_link(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_persistance_length: f64, nondimensional_end_to_end_length_per_link: f64, temperature: f64) -> f64
{
    super::nondimensional_gibbs_free_energy_per_link(&number_of_links, &link_length, &hinge_mass, &nondimensional_persistance_length, &nondimensional_end_to_end_length_per_link, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isometric_nondimensional_relative_gibbs_free_energy(number_of_links: u8, nondimensional_persistance_length: f64, nondimensional_end_to_end_length_per_link: f64) -> f64
{
    super::nondimensional_relative_gibbs_free_energy(&number_of_links, &nondimensional_persistance_length, &nondimensional_end_to_end_length_per_link)
}
#[no_mangle]
pub extern fn physics_single_chain_wlc_thermodynamics_isometric_nondimensional_relative_gibbs_free_energy_per_link(number_of_links: u8, nondimensional_persistance_length: f64, nondimensional_end_to_end_length_per_link: f64) -> f64
{
    super::nondimensional_relative_gibbs_free_energy_per_link(&number_of_links, &nondimensional_persistance_length, &nondimensional_end_to_end_length_per_link)
}