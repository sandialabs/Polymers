#[no_mangle]
pub extern fn fjc_thermodynamics_isometric_force(number_of_links: u8, link_length: f64, hinge_mass: f64, end_to_end_length_per_link: f64, temperature: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).force(&end_to_end_length_per_link, &temperature)
}
#[no_mangle]
pub extern fn fjc_thermodynamics_isometric_nondimensional_force(number_of_links: u8, link_length: f64, hinge_mass: f64, nondimensional_end_to_end_length_per_link: f64) -> f64
{
    super::FJC::init(number_of_links, link_length, hinge_mass).nondimensional_force(&nondimensional_end_to_end_length_per_link)
}