// #[cfg(feature = "extern")]
// pub mod ex;

#[cfg(feature = "python")]
pub mod py;

use crate::physics::single_chain::efjc::thermodynamics::isometric::asymptotic::reduced::legendre::
{
    nondimensional_force,
    nondimensional_relative_helmholtz_free_energy
};
use std::
{
    array::from_fn,
    f64::consts::PI
};

const NUMGRID: usize = 256;

/// The Buche-Silberstein hyperelastic damage constitutive model.
pub struct BucheSilberstein
{
    element: [[f64; NUMGRID]; NUMGRID],
    factor: f64,
    grid: [f64; NUMGRID],
    method: u8,
    normalization: f64,

    /// The nondimensional critical extension which irreversibly breaks chains.
    pub nondimensional_critical_extension: f64,

    /// The nondimensional stiffness of each link in a chain.
    pub nondimensional_link_stiffness: f64,

    /// The number of links in a chain.
    pub number_of_links: u8,

    /// The volumetric swelling ratio.
    pub swelling_ratio: f64
}

pub fn init(method: u8, nondimensional_critical_extension: f64, nondimensional_link_stiffness: f64, number_of_links: u8, swelling_ratio: f64) -> ([[f64; NUMGRID]; NUMGRID], f64, [f64; NUMGRID], f64)
{
    let grid: [f64; NUMGRID] = from_fn(|i|
        nondimensional_critical_extension * (i + 1) as f64 / NUMGRID as f64
    );
    let dz = grid[1] - grid[0];
    let mut element = [[0.0_f64; NUMGRID]; NUMGRID];
    element.iter_mut().zip(grid.iter()).for_each(|(element_i, z_i)|
        element_i.iter_mut().zip(grid.iter()).for_each(|(element_ij, r_j)|
            *element_ij = nondimensional_force(&nondimensional_link_stiffness, &(z_i.powi(2) + r_j.powi(2)).sqrt())
                * (2.0 * z_i.powi(2) - r_j.powi(2)) * 2.0 * PI * (number_of_links as f64) * r_j
                / (z_i.powi(2) + r_j.powi(2)).sqrt() * (dz / swelling_ratio).powi(2)
        )
    );
    let factor;
    let normalization;
    match method {
        1 => panic!("Method number 1 not yet implemented."),
        2 => {
            factor = 0.0;
            normalization = grid.iter().map(|gamma_i|
                gamma_i.powi(2) / nondimensional_relative_helmholtz_free_energy(
                    &number_of_links, &nondimensional_link_stiffness, gamma_i
                ).exp()
            ).sum::<f64>() * 4.0 * PI * dz / swelling_ratio;
        }
        3 => {
            factor = 1.5 * number_of_links as f64
                * nondimensional_link_stiffness * (nondimensional_link_stiffness + 1.0)
                / (nondimensional_link_stiffness.powi(2) + 6.0 * nondimensional_link_stiffness + 3.0);
            normalization = grid.iter().map(|gamma_i|
                gamma_i.powi(2) / (factor * gamma_i.powi(2)).exp()
            ).sum::<f64>() * 4.0 * PI * dz / swelling_ratio;
        }
        _ => panic!("Invalid method number.")
    }
    (element, factor, grid, normalization)
}

#[allow(clippy::too_many_arguments)]
pub fn uniaxial_tension(element: &[[f64; NUMGRID]], factor: &f64, grid: &[f64], normalization: &f64, method: &u8, nondimensional_critical_extension: &f64, nondimensional_link_stiffness: &f64, number_of_links: &u8, swelling_ratio: &f64, stretch: &f64, maximum_previous_stretch: &f64) -> [f64; 2]
{
    assert!(stretch >= &1.0);
    let cauchy_stress: f64;
    let total_probability: f64;
    let j_1_3 = swelling_ratio.powf(1.0/3.0);
    match method {
        1 => panic!("Method number 1 not yet implemented."),
        2 => {
            cauchy_stress = element.iter().zip(grid.iter()).flat_map(|(element_i, z_i)|
                element_i.iter().zip(grid.iter()).map(|(element_ij, r_j)|
                    element_ij / nondimensional_relative_helmholtz_free_energy(
                        number_of_links, nondimensional_link_stiffness, &((*z_i / stretch).powi(2) + stretch * r_j.powi(2)).sqrt()
                    ).exp() * ((
                        ((*z_i * maximum_previous_stretch / stretch).powi(2) + stretch / maximum_previous_stretch * r_j.powi(2)).sqrt()
                        <= *nondimensional_critical_extension
                    ) as u8 as f64)
                )
            ).sum::<f64>() / normalization;
            total_probability = grid.iter().flat_map(|z_i|
                grid.iter().map(|r_j|
                    r_j / nondimensional_relative_helmholtz_free_energy(
                        number_of_links, nondimensional_link_stiffness, &((*z_i / stretch).powi(2) + stretch * r_j.powi(2)).sqrt()
                    ).exp() * ((
                        ((*z_i * maximum_previous_stretch / stretch).powi(2) + stretch / maximum_previous_stretch * r_j.powi(2)).sqrt()
                        <= *nondimensional_critical_extension
                    ) as u8 as f64)
                )
            ).sum::<f64>() / normalization * 4.0 * PI * ((grid[1] - grid[0]) / swelling_ratio).powi(2);
        }
        3 => {
            cauchy_stress = element.iter().zip(grid.iter()).flat_map(|(element_i, z_i)|
                element_i.iter().zip(grid.iter()).map(|(element_ij, r_j)|
                    element_ij / (factor * ((*z_i / stretch).powi(2) + stretch * r_j.powi(2)) / j_1_3.powi(2)).exp() * ((
                        ((*z_i * maximum_previous_stretch / stretch).powi(2) + stretch / maximum_previous_stretch * r_j.powi(2)).sqrt()
                        <= *nondimensional_critical_extension
                    ) as u8 as f64)
                )
            ).sum::<f64>() / normalization;
            total_probability = grid.iter().flat_map(|z_i|
                grid.iter().map(|r_j|
                    r_j / (factor * ((*z_i / stretch).powi(2) + stretch * r_j.powi(2)) / j_1_3.powi(2)).exp() * ((
                        ((*z_i * maximum_previous_stretch / stretch).powi(2) + stretch / maximum_previous_stretch * r_j.powi(2)).sqrt()
                        <= *nondimensional_critical_extension
                    ) as u8 as f64)
                )
            ).sum::<f64>() / normalization * 4.0 * PI * ((grid[1] - grid[0]) / swelling_ratio).powi(2);
        }
        _ => panic!("Invalid method number.")
    }
    [cauchy_stress, total_probability]
}

/// The implemented functionality of the Buche-Silberstein hyperelastic damage constitutive model.
impl BucheSilberstein
{
    /// Initializes and returns an instance of the model.
    /// 1. Helmholtz method for both the Helmholtz free energy and the equilibrium distribution.
    /// 2. Gibbs-Legendre method for both the Helmholtz free energy and the equilibrium distribution.
    /// 3. Gibbs-Legendre for the Helmholtz free energy and a Gaussian equilibrium distribution.
    pub fn init(method: u8, nondimensional_critical_extension: f64, nondimensional_link_stiffness: f64, number_of_links: u8, swelling_ratio: f64) -> Self
    {
        let (element, factor, grid, normalization) = init(method, nondimensional_critical_extension, nondimensional_link_stiffness, number_of_links, swelling_ratio);
        BucheSilberstein
        {
            element,
            factor,
            grid,
            method,
            normalization,
            nondimensional_critical_extension,
            nondimensional_link_stiffness,
            number_of_links,
            swelling_ratio
        }
    }
    /// The nondimensional Cauchy stress as a function of stretch in uniaxial tension.
    pub fn uniaxial_tension<const L: usize>(&self, stretch: &[f64]) -> [[f64; 2]; L]
    {
        let mut maximum_previous_stretch = 1.0;
        from_fn(|i|{
            maximum_previous_stretch =
            if maximum_previous_stretch < stretch[i] {
                stretch[i]
            } else {
                maximum_previous_stretch
            };
            uniaxial_tension(&self.element, &self.factor, &self.grid, &self.normalization, &self.method, &self.nondimensional_critical_extension, &self.nondimensional_link_stiffness, &self.number_of_links, &self.swelling_ratio, &stretch[i], &maximum_previous_stretch)
        })
    }
}