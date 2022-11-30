#![cfg(test)]
use super::*;
use crate::physics::BOLTZMANN_CONSTANT;
pub struct Parameters
{
    pub abs_tol: f64,
    pub rel_tol: f64,
    pub log_log_tol: f64,
    pub log_log_scale: f64,
    pub number_of_loops: u32,
    pub hinge_mass_reference: f64,
    pub hinge_mass_scale: f64,
    pub link_length_reference: f64,
    pub link_length_scale: f64,
    pub number_of_links_minimum: u8,
    pub number_of_links_maximum: u8,
    pub nondimensional_end_to_end_length_per_link_reference: f64,
    pub nondimensional_end_to_end_length_per_link_scale: f64,
    pub nondimensional_potential_distance_reference: f64,
    pub nondimensional_potential_distance_scale: f64,
    pub nondimensional_potential_distance_small: f64,
    pub nondimensional_potential_stiffness_reference: f64,
    pub nondimensional_potential_stiffness_scale: f64,
    pub temperature_reference: f64,
    pub temperature_scale: f64,
}
impl Default for Parameters
{
    fn default() -> Self
    {
        Self
        {
            abs_tol: 1e-8,
            rel_tol: 1e-6,
            log_log_tol: 1e-3,
            log_log_scale: 12e-1,
            number_of_loops: 8,
            hinge_mass_reference: 1e0,
            hinge_mass_scale: 1e0,
            link_length_reference: 1e0,
            link_length_scale: 1e0,
            number_of_links_minimum: 4,
            number_of_links_maximum: 25,
            nondimensional_end_to_end_length_per_link_reference: 5e-1,
            nondimensional_end_to_end_length_per_link_scale: 99e-2,
            nondimensional_potential_distance_reference: 1e0,
            nondimensional_potential_distance_scale: 2e0,
            nondimensional_potential_distance_small: 25e-2,
            nondimensional_potential_stiffness_reference: 5e2,
            nondimensional_potential_stiffness_scale: 1e3,
            temperature_reference: 3e2,
            temperature_scale: 1e2,
        }
    }
}
mod base
{
    use super::*;
    use rand::Rng;
    #[test]
    fn init()
    {
        let parameters = Parameters::default();
        let _ = FJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference);
    }
    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            assert_eq!(number_of_links, FJC::init(number_of_links, parameters.link_length_reference, parameters.hinge_mass_reference).number_of_links);
        }
    }
    #[test]
    fn link_length()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(link_length, FJC::init(parameters.number_of_links_minimum, link_length, parameters.hinge_mass_reference).link_length);
        }
    }
    #[test]
    fn hinge_mass()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(hinge_mass, FJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, hinge_mass).hinge_mass);
        }
    }
    #[test]
    fn number_of_links_and_link_length_and_hinge_mass()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_length = rng.gen::<f64>();
            assert_eq!(link_length, FJC::init(number_of_links, link_length, hinge_mass).link_length);
        }
    }
}
// mod legendre
// {
//     test    isometricLegendre == isotensional
//     and  isotensionalLegendre == isometric
// }
// mod thermodynamic_limit
// {
//     test isometric == isotensional
// }
mod strong_potential_small_distance
{
    use super::*;
    use rand::Rng;
    use crate::math::integrate;
    use crate::physics::single_chain::fjc::thermodynamics::modified_canonical::{
        ZERO,
        POINTS
    };
    #[test]
    fn relative_helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
                let integrand_numerator = |end_to_end_length: f64|
                {
                    (model.modified_canonical.relative_helmholtz_free_energy(&end_to_end_length, &potential_stiffness, &temperature) - model.isometric.relative_helmholtz_free_energy(&end_to_end_length, &temperature)).powf(2.0)
                };
                let integrand_denominator = |end_to_end_length: f64|
                {
                    model.modified_canonical.relative_helmholtz_free_energy(&end_to_end_length, &potential_stiffness, &temperature).powf(2.0)
                };
                let numerator = integrate(integrand_numerator, ZERO, parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length, POINTS);
                let denominator = integrate(integrand_denominator, ZERO, parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length, POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_scale);
            let residual_rel_2 = residual_rel(parameters.log_log_scale*parameters.nondimensional_potential_stiffness_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!((log_log_slope + 1.0).abs()  <= parameters.log_log_tol);
        }
    }
    #[test]
    fn relative_helmholtz_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
                let integrand_numerator = |end_to_end_length: f64|
                {
                    (model.modified_canonical.relative_helmholtz_free_energy_per_link(&end_to_end_length, &potential_stiffness, &temperature) - model.isometric.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature)).powf(2.0)
                };
                let integrand_denominator = |end_to_end_length: f64|
                {
                    model.modified_canonical.relative_helmholtz_free_energy_per_link(&end_to_end_length, &potential_stiffness, &temperature).powf(2.0)
                };
                let numerator = integrate(integrand_numerator, ZERO, parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length, POINTS);
                let denominator = integrate(integrand_denominator, ZERO, parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length, POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_scale);
            let residual_rel_2 = residual_rel(parameters.log_log_scale*parameters.nondimensional_potential_stiffness_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!((log_log_slope + 1.0).abs()  <= parameters.log_log_tol);
        }
    }
    #[test]
    fn nondimensional_relative_helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let integrand_numerator = |nondimensional_end_to_end_length_per_link: f64|
                {
                    (model.modified_canonical.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &nondimensional_potential_stiffness) - model.isometric.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link)).powf(2.0)
                };
                let integrand_denominator = |nondimensional_end_to_end_length_per_link: f64|
                {
                    model.modified_canonical.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &nondimensional_potential_stiffness).powf(2.0)
                };
                let numerator = integrate(integrand_numerator, ZERO, parameters.nondimensional_potential_distance_small, POINTS);
                let denominator = integrate(integrand_denominator, ZERO, parameters.nondimensional_potential_distance_small, POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_scale);
            let residual_rel_2 = residual_rel(parameters.log_log_scale*parameters.nondimensional_potential_stiffness_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!((log_log_slope + 1.0).abs()  <= parameters.log_log_tol);
        }
    }
    #[test]
    fn nondimensional_relative_helmholtz_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let integrand_numerator = |nondimensional_end_to_end_length_per_link: f64|
                {
                    (model.modified_canonical.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &nondimensional_potential_stiffness) - model.isometric.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link)).powf(2.0)
                };
                let integrand_denominator = |nondimensional_end_to_end_length_per_link: f64|
                {
                    model.modified_canonical.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &nondimensional_potential_stiffness).powf(2.0)
                };
                let numerator = integrate(integrand_numerator, ZERO, parameters.nondimensional_potential_distance_small, POINTS);
                let denominator = integrate(integrand_denominator, ZERO, parameters.nondimensional_potential_distance_small, POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_scale);
            let residual_rel_2 = residual_rel(parameters.log_log_scale*parameters.nondimensional_potential_stiffness_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!((log_log_slope + 1.0).abs()  <= parameters.log_log_tol);
        }
    }
}
// mod weak_potential_high_distance
// {
//     test isometric == modifiedCanonical
// }