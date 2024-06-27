#[cfg(feature = "extern")]
pub mod ex;

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

pub const NUMGRID: usize = 256;

/// The Buche-Silberstein hyperelastic constitutive model.
pub struct BucheSilberstein
{
    element: [[f64; NUMGRID]; NUMGRID],
    factor: f64,
    grid: [f64; NUMGRID],
    method: u8,
    normalization: f64,

    /// The nondimensional stiffness of each link in a chain.
    pub nondimensional_link_stiffness: f64,

    /// The number of links in a chain.
    pub number_of_links: u8
}

pub fn init(method: u8, nondimensional_link_stiffness: f64, number_of_links: u8) -> ([[f64; NUMGRID]; NUMGRID], f64, [f64; NUMGRID], f64)
{
    let num_grid = NUMGRID as f64 + 1.0;
    let w: [f64; NUMGRID] = from_fn(|i|
        (i + 1) as f64 / num_grid
    );
    let dw = w[1] - w[0];
    let grid: [f64; NUMGRID] = from_fn(|i|
        w[i].atanh()
    );
    let mut element = [[0.0_f64; NUMGRID]; NUMGRID];
    element.iter_mut().zip(grid.iter().zip(w.iter())).for_each(|(element_i, (z_i, w_i))|
        element_i.iter_mut().zip(grid.iter().zip(w.iter())).for_each(|(element_ij, (r_j, w_j))|
            *element_ij = nondimensional_force(&nondimensional_link_stiffness, &(z_i.powi(2) + r_j.powi(2)).sqrt())
                * (2.0 * z_i.powi(2) - r_j.powi(2)) * 2.0 * PI * (number_of_links as f64) * r_j
                / (z_i.powi(2) + r_j.powi(2)).sqrt() / (1.0 - w_i.powi(2)) / (1.0 - w_j.powi(2)) * dw.powi(2)
        )
    );
    let factor;
    let normalization;
    match method {
        1 => panic!("Method number 1 not yet implemented."),
        2 => {
            factor = 0.0;
            normalization = grid.iter().zip(w.iter()).map(|(gamma_i, w_i)|
                gamma_i.powi(2) / (1.0 - w_i.powi(2)) / nondimensional_relative_helmholtz_free_energy(
                    &number_of_links, &nondimensional_link_stiffness, gamma_i
                ).exp()
            ).sum::<f64>() * 4.0 * PI * dw;
        }
        3 => {
            factor = 1.5 * number_of_links as f64
                * nondimensional_link_stiffness * (nondimensional_link_stiffness + 1.0)
                / (nondimensional_link_stiffness.powi(2) + 6.0 * nondimensional_link_stiffness + 3.0);
            normalization = grid.iter().zip(w.iter()).map(|(gamma_i, w_i)|
                gamma_i.powi(2) / (1.0 - w_i.powi(2)) / (factor * gamma_i.powi(2)).exp()
            ).sum::<f64>() * 4.0 * PI * dw;
        }
        _ => panic!("Invalid method number.")
    }
    (element, factor, grid, normalization)
}

#[allow(clippy::too_many_arguments)]
pub fn uniaxial_tension(element: &[[f64; NUMGRID]], factor: &f64, grid: &[f64], normalization: &f64, method: &u8, nondimensional_link_stiffness: &f64, number_of_links: &u8, stretch: &f64) -> f64
{
    match method {
        1 => panic!("Method number 1 not yet implemented."),
        2 => {
            element.iter().zip(grid.iter()).flat_map(|(element_i, z_i)|
                element_i.iter().zip(grid.iter()).map(|(element_ij, r_j)|
                    element_ij / nondimensional_relative_helmholtz_free_energy(
                        number_of_links, nondimensional_link_stiffness, &((*z_i / stretch).powi(2) + stretch * r_j.powi(2)).sqrt()
                    ).exp()
                )
            ).sum::<f64>() / normalization
        }
        3 => {
            element.iter().zip(grid.iter()).flat_map(|(element_i, z_i)|
                element_i.iter().zip(grid.iter()).map(|(element_ij, r_j)|
                    element_ij / (factor * ((*z_i / stretch).powi(2) + stretch * r_j.powi(2))).exp()
                )
            ).sum::<f64>() / normalization
        }
        _ => panic!("Invalid method number.")
    }
}

#[allow(clippy::too_many_arguments)]
pub fn equibiaxial_tension(element: &[[f64; NUMGRID]], factor: &f64, grid: &[f64], normalization: &f64, method: &u8, nondimensional_link_stiffness: &f64, number_of_links: &u8, stretch: &f64) -> f64
{
    match method {
        1 => panic!("Method number 1 not yet implemented."),
        2 => {
            element.iter().zip(grid.iter()).flat_map(|(element_i, z_i)|
                element_i.iter().zip(grid.iter()).map(|(element_ij, r_j)|
                    -element_ij / nondimensional_relative_helmholtz_free_energy(
                        number_of_links, nondimensional_link_stiffness, &((*z_i * stretch.powi(2)).powi(2) + (*r_j / stretch).powi(2)).sqrt()
                    ).exp()
                )
            ).sum::<f64>() / normalization
        }
        3 => {
            element.iter().zip(grid.iter()).flat_map(|(element_i, z_i)|
                element_i.iter().zip(grid.iter()).map(|(element_ij, r_j)|
                    -element_ij / (factor * ((*z_i * stretch.powi(2)).powi(2) + (*r_j / stretch).powi(2))).exp()
                )
            ).sum::<f64>() / normalization
        }
        _ => panic!("Invalid method number.")
    }
}

/// The implemented functionality of the Buche-Silberstein hyperelastic constitutive model.
impl BucheSilberstein
{
    /// Initializes and returns an instance of the model.
    /// 1. Helmholtz method for both the Helmholtz free energy and the equilibrium distribution.
    /// 2. Gibbs-Legendre method for both the Helmholtz free energy and the equilibrium distribution.
    /// 3. Gibbs-Legendre for the Helmholtz free energy and a Gaussian equilibrium distribution.
    pub fn init(method: u8, nondimensional_link_stiffness: f64, number_of_links: u8) -> Self
    {
        let (element, factor, grid, normalization) = init(method, nondimensional_link_stiffness, number_of_links);
        BucheSilberstein
        {
            element,
            factor,
            grid,
            method,
            normalization,
            nondimensional_link_stiffness,
            number_of_links
        }
    }
    /// The nondimensional Cauchy stress as a function of stretch in uniaxial tension.
    pub fn uniaxial_tension(&self, stretch: &f64) -> f64
    {
        uniaxial_tension(&self.element, &self.factor, &self.grid, &self.normalization, &self.method, &self.nondimensional_link_stiffness, &self.number_of_links, stretch)
    }
    /// The nondimensional Cauchy stress as a function of stretch in equibiaxial tension.
    pub fn equibiaxial_tension(&self, stretch: &f64) -> f64
    {
        equibiaxial_tension(&self.element, &self.factor, &self.grid, &self.normalization, &self.method, &self.nondimensional_link_stiffness, &self.number_of_links, stretch)
    }
}