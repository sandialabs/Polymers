mod test;

use rand::prelude::*;

use rand_distr::{Normal, Distribution};

use std::f64::consts::TAU as TWO_PI;

pub fn random_configuration<const NUMBER_OF_LINKS: usize>(dist: Normal<f64>, rng: &mut ThreadRng) -> [[f64; 3]; NUMBER_OF_LINKS]
{
    let mut lambda: f64 = 0.0;
    let mut phi: f64 = 0.0;
    let mut theta: f64 = 0.0;
    let mut position = [0.0; 3];
    let mut configuration = [[0.0; 3]; NUMBER_OF_LINKS];
    configuration.iter_mut().for_each(|coordinate|{
        lambda = dist.sample(rng);
        phi = TWO_PI * rng.gen::<f64>();
        theta = (1.0 - 2.0 * rng.gen::<f64>()).acos();
        position[0] += lambda * (theta.sin() * phi.cos());
        position[1] += lambda * (theta.sin() * phi.sin());
        position[2] += lambda * theta.cos();
        coordinate.iter_mut().zip(position.iter()).for_each(|(coordinate_i, position_i)|
            *coordinate_i = *position_i
        );
    });
    configuration
}

pub fn random_nondimensional_end_to_end_length<const NUMBER_OF_LINKS: usize>(dist: Normal<f64>, rng: &mut ThreadRng) -> f64
{
    random_configuration::<NUMBER_OF_LINKS>(dist, rng)[NUMBER_OF_LINKS - 1].iter().map(|entry| entry * entry).sum::<f64>().sqrt()
}

pub fn nondimensional_equilibrium_radial_distribution<const NUMBER_OF_BINS: usize, const NUMBER_OF_LINKS: usize>(gamma_max: &f64, kappa: &f64, number_of_samples: usize) -> ([f64; NUMBER_OF_BINS], [f64; NUMBER_OF_BINS])
{
    let mut rng = rand::thread_rng();
    let dist = Normal::new(1.0, 1.0/kappa.sqrt()).unwrap();
    let number_of_links_f64 = NUMBER_OF_LINKS as f64;
    let mut bin_edges = [0.0_f64; NUMBER_OF_BINS];
    bin_edges.iter_mut().enumerate().for_each(|(bin_index, bin_edge)|
        *bin_edge = gamma_max * (bin_index as f64 + 1.0)/(NUMBER_OF_BINS as f64)
    );
    let mut bin_counts = [0_u128; NUMBER_OF_BINS];
    let mut gamma: f64 = 0.0;
    (0..number_of_samples).for_each(|_|{
        gamma = random_nondimensional_end_to_end_length::<NUMBER_OF_LINKS>(dist, &mut rng)/number_of_links_f64;
        for (bin_edge, bin_count) in bin_edges.iter().zip(bin_counts.iter_mut())
        {
            if &gamma > gamma_max
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
