#![cfg(test)]

use super::*;

const NUMBER_OF_BINS: usize = 100;
const NUMBER_OF_LINKS: usize = 8;
const NUMBER_OF_SAMPLES: usize = 100000;

#[test]
fn monte_carlo_random_configuration()
{
    let mut rng = rand::thread_rng();
    let configuration = random_configuration::<NUMBER_OF_LINKS>(&mut rng);
    assert_eq!(configuration.len(), NUMBER_OF_LINKS);
    assert_eq!(configuration[0].len(), 3);
}

#[test]
fn monte_carlo_random_nondimensional_end_to_end_length()
{
    let mut rng = rand::thread_rng();
    let gamma = random_nondimensional_end_to_end_length::<NUMBER_OF_LINKS>(&mut rng)/(NUMBER_OF_LINKS as f64);
    assert!(gamma >= 0.0);
    assert!(gamma <= 1.0);
}

#[test]
fn monte_carlo_nondimensional_equilibrium_radial_distribution()
{
    let (gamma, g_eq) = nondimensional_equilibrium_radial_distribution::<NUMBER_OF_BINS, NUMBER_OF_LINKS>(NUMBER_OF_SAMPLES);
    assert_eq!(gamma.len(), NUMBER_OF_BINS);
    assert_eq!(g_eq.len(), NUMBER_OF_BINS);
    gamma.iter().zip(g_eq.iter()).for_each(|(gamma_i, g_eq_i)|{
        assert!(gamma_i > &0.0);
        assert!(gamma_i < &1.0);
        assert!(g_eq_i >= &0.0);
    });
}