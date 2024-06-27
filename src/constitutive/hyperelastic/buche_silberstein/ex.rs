#[no_mangle]
pub unsafe extern fn constitutive_hyperelastic_buche_silberstein_uniaxial_tension(raw_element: *const [[f64; super::NUMGRID]; super::NUMGRID], factor: f64, raw_grid: *const [f64; super::NUMGRID], normalization: f64, method: u8, nondimensional_link_stiffness: f64, number_of_links: u8, stretch: f64) -> f64
{
    let grid = std::slice::from_raw_parts(raw_grid, super::NUMGRID)[0];
    let element = std::slice::from_raw_parts(raw_element, super::NUMGRID * super::NUMGRID)[0];
    super::uniaxial_tension(&element, &factor, &grid, &normalization, &method, &nondimensional_link_stiffness, &number_of_links, &stretch)
}
#[no_mangle]
pub unsafe extern fn constitutive_hyperelastic_buche_silberstein_equibiaxial_tension(raw_element: *const [[f64; super::NUMGRID]; super::NUMGRID], factor: f64, raw_grid: *const [f64; super::NUMGRID], normalization: f64, method: u8, nondimensional_link_stiffness: f64, number_of_links: u8, stretch: f64) -> f64
{
    let grid = std::slice::from_raw_parts(raw_grid, super::NUMGRID)[0];
    let element = std::slice::from_raw_parts(raw_element, super::NUMGRID * super::NUMGRID)[0];
    super::equibiaxial_tension(&element, &factor, &grid, &normalization, &method, &nondimensional_link_stiffness, &number_of_links, &stretch)
}