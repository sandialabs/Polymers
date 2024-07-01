mod test;

use rand::prelude::*;

use std::f64::consts::PI;
use std::f64::consts::TAU as TWO_PI;

pub fn cross(u: &[f64; 3], v: &[f64; 3]) -> [f64; 3]
{
    [u[1] * v[2] - u[2] * v[1],
     u[2] * v[0] - u[0] * v[2],
     u[0] * v[1] - u[1] * v[0]]
}

pub fn random_configuration<const NUMBER_OF_LINKS: usize>(theta: &f64, rng: &mut ThreadRng) -> [[f64; 3]; NUMBER_OF_LINKS]
{
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
    configuration.iter_mut().for_each(|coordinate|{
        phi = TWO_PI * rng.gen::<f64>();
        phi_cos = phi.cos();
        phi_sin = phi.sin();
        t = std::array::from_fn(|_| rng.gen::<f64>());
        u = cross(&r, &t);
        u_normalization = u.iter().map(|u_i| u_i * u_i).sum::<f64>().sqrt();
        u.iter_mut().for_each(|u_i| *u_i /= u_normalization);
        v = cross(&r, &u);
        r.iter_mut().zip(u.iter().zip(v.iter())).for_each(|(r_i, (u_i, v_i))|
            *r_i = (u_i * phi_cos + v_i * phi_sin) * theta_sin + *r_i * theta_cos
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

pub fn random_nondimensional_end_to_end_length<const NUMBER_OF_LINKS: usize>(theta: &f64, rng: &mut ThreadRng) -> f64
{
    random_configuration::<NUMBER_OF_LINKS>(theta, rng)[NUMBER_OF_LINKS - 1].iter().map(|entry| entry * entry).sum::<f64>().sqrt()
}

pub fn nondimensional_equilibrium_radial_distribution<const NUMBER_OF_BINS: usize, const NUMBER_OF_LINKS: usize>(theta: &f64, number_of_samples: usize) -> ([f64; NUMBER_OF_BINS], [f64; NUMBER_OF_BINS])
{
    let mut rng = rand::thread_rng();
    let number_of_links_f64 = NUMBER_OF_LINKS as f64;
    let gamma_max = (2.0 - 2.0*(PI - theta).cos()).sqrt()/2.0;
    let mut bin_edges = [0.0_f64; NUMBER_OF_BINS];
    bin_edges.iter_mut().enumerate().for_each(|(bin_index, bin_edge)|
        *bin_edge = gamma_max * (bin_index as f64 + 1.0)/(NUMBER_OF_BINS as f64)
    );
    let mut bin_counts = [0_u128; NUMBER_OF_BINS];
    let mut gamma: f64 = 0.0;
    (0..number_of_samples).for_each(|_|{
        gamma = random_nondimensional_end_to_end_length::<NUMBER_OF_LINKS>(theta, &mut rng)/number_of_links_f64;
        for (bin_edge, bin_count) in bin_edges.iter().zip(bin_counts.iter_mut())
        {
            if gamma > gamma_max
            {
                panic!()
            }
            if &gamma < bin_edge
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
    let mut bin_centers = [0.0_f64; NUMBER_OF_BINS];
    bin_centers.iter_mut().enumerate().for_each(|(bin_index, bin_center)|
        *bin_center = gamma_max * (bin_index as f64 + 0.5)/(NUMBER_OF_BINS as f64)
    );
    (bin_centers, bin_probabilities)
}