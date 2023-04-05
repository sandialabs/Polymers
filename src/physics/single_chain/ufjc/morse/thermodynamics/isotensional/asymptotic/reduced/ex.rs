#[no_mangle]
pub extern fn physics_single_chain_ufjc_morse_thermodynamics_isotensional_asymptotic_reduced_end_to_end_length(number_of_links: u8, link_length: f64, link_stiffness: f64, link_energy: f64, force: f64, temperature: f64) -> f64
{
    super::end_to_end_length(&number_of_links, &link_length, &link_stiffness, &link_energy, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_morse_thermodynamics_isotensional_asymptotic_reduced_end_to_end_length_per_link(link_length: f64, link_stiffness: f64, link_energy: f64, force: f64, temperature: f64) -> f64
{
    super::end_to_end_length_per_link(&link_length, &link_stiffness, &link_energy, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_morse_thermodynamics_isotensional_asymptotic_reduced_nondimensional_end_to_end_length(number_of_links: u8, nondimensional_link_stiffness: f64, nondimensional_link_energy: f64, nondimensional_force: f64) -> f64
{
    super::nondimensional_end_to_end_length(&number_of_links, &nondimensional_link_stiffness, &nondimensional_link_energy, &nondimensional_force)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_morse_thermodynamics_isotensional_asymptotic_reduced_nondimensional_end_to_end_length_per_link(nondimensional_link_stiffness: f64, nondimensional_link_energy: f64, nondimensional_force: f64) -> f64
{
    super::nondimensional_end_to_end_length_per_link(&nondimensional_link_stiffness, &nondimensional_link_energy, &nondimensional_force)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_morse_thermodynamics_isotensional_asymptotic_reduced_gibbs_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, link_stiffness: f64, link_energy: f64, force: f64, temperature: f64) -> f64
{
    super::gibbs_free_energy(&number_of_links, &link_length, &hinge_mass, &link_stiffness, &link_energy, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_morse_thermodynamics_isotensional_asymptotic_reduced_gibbs_free_energy_per_link(link_length: f64, hinge_mass: f64, link_stiffness: f64, link_energy: f64, force: f64, temperature: f64) -> f64
{
    super::gibbs_free_energy_per_link(&link_length, &hinge_mass, &link_stiffness, &link_energy, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_morse_thermodynamics_isotensional_asymptotic_reduced_relative_gibbs_free_energy(number_of_links: u8, link_length: f64, link_stiffness: f64, link_energy: f64, force: f64, temperature: f64) -> f64
{
    super::relative_gibbs_free_energy(&number_of_links, &link_length, &link_stiffness, &link_energy, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_morse_thermodynamics_isotensional_asymptotic_reduced_relative_gibbs_free_energy_per_link(link_length: f64, link_stiffness: f64, link_energy: f64, force: f64, temperature: f64) -> f64
{
    super::relative_gibbs_free_energy_per_link(&link_length, &link_stiffness, &link_energy, &force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_morse_thermodynamics_isotensional_asymptotic_reduced_nondimensional_gibbs_free_energy(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_link_stiffness: f64, nondimensional_link_energy: f64, nondimensional_force: f64, temperature: f64) -> f64
{
    super::nondimensional_gibbs_free_energy(&number_of_links, &link_length, &hinge_mass, &nondimensional_link_stiffness, &nondimensional_link_energy, &nondimensional_force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_morse_thermodynamics_isotensional_asymptotic_reduced_nondimensional_gibbs_free_energy_per_link(link_length: f64, hinge_mass: f64, nondimensional_link_stiffness: f64, nondimensional_link_energy: f64, nondimensional_force: f64, temperature: f64) -> f64
{
    super::nondimensional_gibbs_free_energy_per_link(&link_length, &hinge_mass, &nondimensional_link_stiffness, &nondimensional_link_energy, &nondimensional_force, &temperature)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_morse_thermodynamics_isotensional_asymptotic_reduced_nondimensional_relative_gibbs_free_energy(number_of_links: u8, nondimensional_link_stiffness: f64, nondimensional_link_energy: f64, nondimensional_force: f64) -> f64
{
    super::nondimensional_relative_gibbs_free_energy(&number_of_links, &nondimensional_link_stiffness, &nondimensional_link_energy, &nondimensional_force)
}
#[no_mangle]
pub extern fn physics_single_chain_ufjc_morse_thermodynamics_isotensional_asymptotic_reduced_nondimensional_relative_gibbs_free_energy_per_link(nondimensional_link_stiffness: f64, nondimensional_link_energy: f64, nondimensional_force: f64) -> f64
{
    super::nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_link_stiffness, &nondimensional_link_energy, &nondimensional_force)
}
