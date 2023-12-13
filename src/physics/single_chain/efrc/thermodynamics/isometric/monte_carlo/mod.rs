mod test;

use rand::prelude::*;

use rand_distr::{Normal, Distribution};

use std::f64::consts::TAU as TWO_PI;

use crate::physics::single_chain::frc::thermodynamics::isometric::monte_carlo::cross;

pub fn random_configuration<const NUMBER_OF_LINKS: usize>(theta: &f64, dist: Normal<f64>, rng: &mut ThreadRng) -> [[f64; 3]; NUMBER_OF_LINKS]
{
    let mut lambda: f64 = 0.0;
    let mut phi: f64 = 0.0;
    let mut phi_cos: f64 = 0.0;
    let mut phi_sin: f64 = 0.0;
    let theta_cos: f64 = theta.cos();
    let theta_sin: f64 = theta.sin();
    let mut configuration = [[0.0; 3]; NUMBER_OF_LINKS];
    let mut position = [0.0; 3];
    let mut r = [1.0, 0.0, 0.0];
    let mut t = [0.0; 3];
    let mut u = [0.0; 3];
    let mut v = [0.0; 3];
    let mut u_normalization = 0.0;
    let mut r_n = [1.0, 0.0, 0.0];
    let mut r_normalization = 0.0;
    configuration.iter_mut().for_each(|coordinate|{
        lambda = dist.sample(rng);
        r_n = r;
        r_normalization = r_n.iter().map(|r_n_i| r_n_i * r_n_i).sum::<f64>().sqrt();
        r_n.iter_mut().for_each(|r_n_i| *r_n_i /= r_normalization);
        phi = TWO_PI * rng.gen::<f64>();
        phi_cos = phi.cos();
        phi_sin = phi.sin();
        t = std::array::from_fn(|_| rng.gen::<f64>());
        u = cross(&r_n, &t);
        u_normalization = u.iter().map(|u_i| u_i * u_i).sum::<f64>().sqrt();
        u.iter_mut().for_each(|u_i| *u_i /= u_normalization);
        v = cross(&r_n, &u);
        r.iter_mut().zip(r_n.iter().zip(u.iter().zip(v.iter()))).for_each(|(r_i, (r_n_i, (u_i, v_i)))|
            *r_i = lambda * ((u_i * phi_cos + v_i * phi_sin) * theta_sin + *r_n_i * theta_cos)
        );
        position.iter_mut().zip(r.iter()).for_each(|(position_i, r_i)|
            *position_i += r_i
        );
        coordinate.iter_mut().zip(position.iter()).for_each(|(coordinate_i, position_i)|
            *coordinate_i = *position_i
        );
    });
    configuration
}

pub fn random_nondimensional_end_to_end_length<const NUMBER_OF_LINKS: usize>(theta: &f64, dist: Normal<f64>, rng: &mut ThreadRng) -> f64
{
    random_configuration::<NUMBER_OF_LINKS>(theta, dist, rng)[NUMBER_OF_LINKS - 1].iter().map(|entry| entry * entry).sum::<f64>().sqrt()
}

pub fn nondimensional_equilibrium_radial_distribution<const NUMBER_OF_BINS: usize, const NUMBER_OF_LINKS: usize>(gamma_max: &f64, kappa: &f64, theta: &f64, number_of_samples: usize) -> ([f64; NUMBER_OF_BINS], [f64; NUMBER_OF_BINS])
{
    let mut rng = rand::thread_rng();
    let dist = Normal::new(1.0, 1.0/kappa.sqrt()).unwrap();
    let number_of_links_f64 = NUMBER_OF_LINKS as f64;
    let mut bin_centers = [0.0_f64; NUMBER_OF_BINS];
    bin_centers.iter_mut().enumerate().for_each(|(bin_index, bin_center)|
        *bin_center = gamma_max * ((bin_index as f64) + 0.5)/(NUMBER_OF_BINS as f64)
    );
    let mut bin_counts = [0_u128; NUMBER_OF_BINS];
    let mut gamma: f64 = 0.0;
    (0..number_of_samples).for_each(|_|{
        gamma = random_nondimensional_end_to_end_length::<NUMBER_OF_LINKS>(theta, dist, &mut rng)/number_of_links_f64;
        for (bin_center, bin_count) in bin_centers.iter().zip(bin_counts.iter_mut())
        {
            if &gamma < bin_center
            {
                *bin_count += 1;
                break
            }
        }
    });
    let normalization = gamma_max * (number_of_samples as f64)/(NUMBER_OF_BINS as f64);
    let mut bin_probabilities = [0.0_f64; NUMBER_OF_BINS];
    bin_probabilities.iter_mut().zip(bin_counts.iter()).for_each(|(bin_probability, bin_count)|
        *bin_probability = (*bin_count as f64)/normalization
    );
    (bin_centers, bin_probabilities)
}