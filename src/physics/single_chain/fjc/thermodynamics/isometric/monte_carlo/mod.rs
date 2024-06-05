mod test;

use rand::prelude::*;

use std::f64::consts::TAU as TWO_PI;

pub fn random_configuration<const NUMBER_OF_LINKS: usize>(rng: &mut ThreadRng) -> ([[f64; 3]; NUMBER_OF_LINKS], [f64; NUMBER_OF_LINKS])
{
    let mut phi: f64 = 0.0;
    let mut theta: f64 = 0.0;
    let mut position = [0.0; 3];
    let mut configuration = [[0.0; 3]; NUMBER_OF_LINKS];
    let mut cos_angles = [0.0; NUMBER_OF_LINKS];
    configuration.iter_mut().zip(cos_angles.iter_mut()).for_each(|(coordinate, cos_angle)|{
        phi = TWO_PI * rng.gen::<f64>();
        theta = (1.0 - 2.0 * rng.gen::<f64>()).acos();
        position[0] += theta.sin() * phi.cos();
        position[1] += theta.sin() * phi.sin();
        position[2] += theta.cos();
        *cos_angle = theta.cos();
        coordinate.iter_mut().zip(position.iter()).for_each(|(coordinate_i, position_i)|
            *coordinate_i = *position_i
        );
    });
    (configuration, cos_angles)
}

pub fn random_nondimensional_end_to_end_length<const NUMBER_OF_LINKS: usize>(rng: &mut ThreadRng) -> (f64, [f64; NUMBER_OF_LINKS])
{
    let (configuration, cos_angles) = random_configuration::<NUMBER_OF_LINKS>(rng);
    let xi = configuration[NUMBER_OF_LINKS - 1].iter().map(|entry| entry * entry).sum::<f64>().sqrt();
    (xi, cos_angles)
}

pub fn nondimensional_equilibrium_radial_distribution<const NUMBER_OF_BINS: usize, const NUMBER_OF_LINKS: usize>(number_of_samples: usize) -> ([f64; NUMBER_OF_BINS], [f64; NUMBER_OF_BINS], [[f64; NUMBER_OF_LINKS]; NUMBER_OF_BINS], [[f64; NUMBER_OF_LINKS]; NUMBER_OF_BINS], [[f64; NUMBER_OF_LINKS]; NUMBER_OF_BINS])
{
    //
    // need to calculate running average of 3 moments for each link in each bin
    //
    // then make test where you integrate for isotensional ensemble and compare to exact
    // (use gamma_0 from exact FJC and corrections from this)
    //
    let mut rng = rand::thread_rng();
    let number_of_links_f64 = NUMBER_OF_LINKS as f64;
    let mut bin_edges = [0.0_f64; NUMBER_OF_BINS];
    bin_edges.iter_mut().enumerate().for_each(|(bin_index, bin_edge)|
        *bin_edge = (bin_index as f64 + 1.0)/(NUMBER_OF_BINS as f64)
    );
    let mut germa: f64 = 0.0;
    let mut gamma: f64 = 0.0;
    let mut cos_angles = [0.0_f64; NUMBER_OF_LINKS];
    // let first_moment: f64 = 0.0;
    let mut bin_counts = [0_u128; NUMBER_OF_BINS];
    let mut first_moments = [[0.0; NUMBER_OF_LINKS]; NUMBER_OF_BINS];
    let second_moments = [[0.0; NUMBER_OF_LINKS]; NUMBER_OF_BINS];
    let third_moments = [[0.0; NUMBER_OF_LINKS]; NUMBER_OF_BINS];
    (0..number_of_samples).for_each(|_sample_num|{
        (germa, cos_angles) = random_nondimensional_end_to_end_length::<NUMBER_OF_LINKS>(&mut rng);
        gamma = germa/number_of_links_f64;
        for (bin_edge, (bin_count, first_moment_links)) in bin_edges.iter().zip(bin_counts.iter_mut().zip(first_moments.iter_mut()))
        {
            if &gamma < bin_edge
            {
                *bin_count += 1;
                first_moment_links.iter_mut().zip(cos_angles.iter()).for_each(|(first_moment, cos_angle)|
                    *first_moment += (cos_angle - *first_moment)/(*bin_count as f64)
                );
                break
            }
        }
    });
    let normalization = (number_of_samples as f64)/(NUMBER_OF_BINS as f64);
    let mut bin_probabilities = [0.0_f64; NUMBER_OF_BINS];
    bin_probabilities.iter_mut().zip(bin_counts.iter()).for_each(|(bin_probability, bin_count)|
        *bin_probability = (*bin_count as f64)/normalization
    );
    let mut bin_centers = [0.0_f64; NUMBER_OF_BINS];
    bin_centers.iter_mut().enumerate().for_each(|(bin_index, bin_center)|
        *bin_center = (bin_index as f64 + 0.5)/(NUMBER_OF_BINS as f64)
    );
    (bin_centers, bin_probabilities, first_moments, second_moments, third_moments)
}
