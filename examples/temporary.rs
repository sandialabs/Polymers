use polymers::physics::single_chain::efjc::thermodynamics::isometric::asymptotic::reduced::legendre::{
    nondimensional_force,
    // nondimensional_relative_helmholtz_free_energy
};
use std::{
    array::from_fn,
    f64::consts::PI
};

const LINKS: u8 = 25;
const KAPPA: f64 = 50.0;
const NUMGRID: usize = 512;

const FAC: f64 = 1.5 * LINKS as f64 * KAPPA * (KAPPA + 1.0) / (KAPPA * KAPPA + 6.0 * KAPPA + 3.0);

fn main() {
    let num_grid = NUMGRID as f64 + 1.0;
    let w: [f64; NUMGRID] = from_fn(|i|
        (i + 1) as f64 / num_grid
    );
    let dw = w[1] - w[0];
    let z: [f64; NUMGRID] = from_fn(|i|
        w[i].atanh()
    );
    let mut element = [[0.0_f64; NUMGRID]; NUMGRID];
    element.iter_mut().zip(z.iter().zip(w.iter())).for_each(|(element_i, (z_i, w_i))|
        element_i.iter_mut().zip(z.iter().zip(w.iter())).for_each(|(element_ij, (r_j, w_j))|
            *element_ij = nondimensional_force(&KAPPA, &(z_i.powi(2) + r_j.powi(2)).sqrt()) * (2.0 * z_i.powi(2) - r_j.powi(2)) * 2.0 * PI * (LINKS as f64) * r_j / (z_i.powi(2) + r_j.powi(2)).sqrt() / (1.0 - w_i.powi(2)) / (1.0 - w_j.powi(2)) * dw.powi(2)
        )
    );
    let normalization: f64 = z.iter().zip(w.iter()).map(|(gamma_i, w_i)|
        // gamma_i.powi(2) / (1.0 - w_i.powi(2)) / nondimensional_relative_helmholtz_free_energy(&LINKS, &KAPPA, &gamma_i).exp()
        gamma_i.powi(2) / (1.0 - w_i.powi(2)) / (FAC * gamma_i.powi(2)).exp()
    ).sum::<f64>() * 4.0 * PI * dw;
    let mut lambda_11: f64;
    let mut sigma_11: f64;
    let num: usize = 33;
    for index in 0..num {
        lambda_11 = 1.0 + 2.5 * (index as f64) / (num as f64 - 1.0);
        sigma_11 = element.iter().zip(z.iter()).flat_map(|(element_i, z_i)|
            element_i.iter().zip(z.iter()).map(|(element_ij, r_j)|
                // element_ij / nondimensional_relative_helmholtz_free_energy(&LINKS, &KAPPA, &((*z_i / lambda_11).powi(2) + lambda_11 * r_j.powi(2)).sqrt()).exp()
                element_ij / (FAC * ((*z_i / lambda_11).powi(2) + lambda_11 * r_j.powi(2))).exp()
            )
        ).sum::<f64>() / normalization;
        println!("{:?}", (lambda_11, sigma_11))
    }
    // gauss is way faster if can justify using it, make options either way
    // seems like GL-GL reproduces past results, but seems like GL-Gauss is still off, why?
}