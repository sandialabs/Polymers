#![cfg(test)]
use super::*;
use crate::physics::single_chain::test::Parameters as DefaultParameters;
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
    pub link_stiffness_reference: f64,
    pub link_stiffness_scale: f64,
    pub nondimensional_link_stiffness_large: f64,
    pub nondimensional_link_stiffness_medium: f64,
    pub nondimensional_force_reference: f64,
    pub nondimensional_force_scale: f64,
    pub temperature_reference: f64,
    pub temperature_scale: f64,
}
impl Default for Parameters
{
    fn default() -> Self
    {
        Self
        {
            number_of_loops: 888,
            abs_tol: DefaultParameters::default().abs_tol,
            rel_tol: DefaultParameters::default().rel_tol,
            log_log_tol: DefaultParameters::default().log_log_tol,
            log_log_scale: DefaultParameters::default().log_log_scale,
            hinge_mass_reference: DefaultParameters::default().hinge_mass_reference,
            hinge_mass_scale: DefaultParameters::default().hinge_mass_scale,
            link_length_reference: DefaultParameters::default().link_length_reference,
            link_length_scale: DefaultParameters::default().link_length_scale,
            number_of_links_minimum: DefaultParameters::default().number_of_links_minimum,
            number_of_links_maximum: DefaultParameters::default().number_of_links_maximum,
            link_stiffness_reference: DefaultParameters::default().link_stiffness_reference,
            link_stiffness_scale: DefaultParameters::default().link_stiffness_scale,
            nondimensional_link_stiffness_large: DefaultParameters::default().nondimensional_link_stiffness_large,
            nondimensional_link_stiffness_medium: DefaultParameters::default().nondimensional_link_stiffness_medium,
            nondimensional_force_reference: DefaultParameters::default().nondimensional_force_reference,
            nondimensional_force_scale: DefaultParameters::default().nondimensional_force_scale,
            temperature_reference: DefaultParameters::default().temperature_reference,
            temperature_scale: DefaultParameters::default().temperature_scale,
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
        let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference);
    }
    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            assert_eq!(number_of_links, EFJC::init(number_of_links, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).number_of_links);
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
            assert_eq!(link_length, EFJC::init(parameters.number_of_links_minimum, link_length, parameters.hinge_mass_reference, parameters.link_stiffness_reference).link_length);
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
            assert_eq!(hinge_mass, EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, hinge_mass, parameters.link_stiffness_reference).hinge_mass);
        }
    }
    #[test]
    fn link_stiffness()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(link_stiffness, EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, link_stiffness).link_stiffness);
        }
    }
    #[test]
    fn number_of_links_and_link_length_and_hinge_mass_and_link_stiffness()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_length = rng.gen::<f64>();
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(link_length, EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness).link_length);
        }
    }
}
mod nondimensional
{
    use super::*;
    use rand::Rng;
    #[test]
    fn end_to_end_length()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_force, &temperature);
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.end_to_end_length(&force, &temperature);
            let residual_abs = &end_to_end_length/link_length - &nondimensional_end_to_end_length;
            let residual_rel = &residual_abs/&nondimensional_end_to_end_length;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn end_to_end_length_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature);
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length_per_link = model.end_to_end_length_per_link(&force, &temperature);
            let residual_abs = &end_to_end_length_per_link/link_length - &nondimensional_end_to_end_length_per_link;
            let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature);
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let gibbs_free_energy = model.gibbs_free_energy(&force, &temperature);
            let residual_abs = &gibbs_free_energy/BOLTZMANN_CONSTANT/temperature - &nondimensional_gibbs_free_energy;
            let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature);
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&force, &temperature);
            let residual_abs = &gibbs_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn relative_gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force, &temperature);
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&force, &temperature);
            let residual_abs = &relative_gibbs_free_energy/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_gibbs_free_energy;
            let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn relative_gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force, &temperature);
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&force, &temperature);
            let residual_abs = &relative_gibbs_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
}
mod per_link
{
    use super::*;
    use rand::Rng;
    #[test]
    fn end_to_end_length()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.end_to_end_length(&force, &temperature);
            let end_to_end_length_per_link = model.end_to_end_length_per_link(&force, &temperature);
            let residual_abs = &end_to_end_length/(number_of_links as f64) - &end_to_end_length_per_link;
            let residual_rel = &residual_abs/&end_to_end_length_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn nondimensional_end_to_end_length()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_force, &temperature);
            let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature);
            let residual_abs = &nondimensional_end_to_end_length/(number_of_links as f64) - &nondimensional_end_to_end_length_per_link;
            let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn gibbs_free_energy()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let gibbs_free_energy = model.gibbs_free_energy(&force, &temperature);
            let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&force, &temperature);
            let residual_abs = &gibbs_free_energy/(number_of_links as f64) - &gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&gibbs_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn relative_gibbs_free_energy()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&force, &temperature);
            let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&force, &temperature);
            let residual_abs = &relative_gibbs_free_energy/(number_of_links as f64) - &relative_gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&relative_gibbs_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn nondimensional_gibbs_free_energy()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature);
            let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature);
            let residual_abs = &nondimensional_gibbs_free_energy/(number_of_links as f64) - &nondimensional_gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn nondimensional_relative_gibbs_free_energy()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force, &temperature);
            let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force, &temperature);
            let residual_abs = &nondimensional_relative_gibbs_free_energy/(number_of_links as f64) - &nondimensional_relative_gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
}
mod relative
{
    use super::*;
    use rand::Rng;
    #[test]
    fn gibbs_free_energy()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let gibbs_free_energy = model.gibbs_free_energy(&force, &temperature);
            let gibbs_free_energy_0 = model.gibbs_free_energy(&(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), &temperature);
            let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&force, &temperature);
            let residual_abs = &gibbs_free_energy - &gibbs_free_energy_0 - &relative_gibbs_free_energy;
            let residual_rel = &residual_abs/&gibbs_free_energy_0;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn gibbs_free_energy_per_link()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&force, &temperature);
            let gibbs_free_energy_per_link_0 = model.gibbs_free_energy_per_link(&(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), &temperature);
            let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&force, &temperature);
            let residual_abs = &gibbs_free_energy_per_link - &gibbs_free_energy_per_link_0 - &relative_gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&gibbs_free_energy_per_link_0;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn nondimensional_gibbs_free_energy()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature);
            let nondimensional_gibbs_free_energy_0 = model.nondimensional_gibbs_free_energy(&ZERO, &temperature);
            let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force, &temperature);
            let residual_abs = &nondimensional_gibbs_free_energy - &nondimensional_gibbs_free_energy_0 - &nondimensional_relative_gibbs_free_energy;
            let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_0;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn nondimensional_gibbs_free_energy_per_link()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature);
            let nondimensional_gibbs_free_energy_per_link_0 = model.nondimensional_gibbs_free_energy_per_link(&ZERO, &temperature);
            let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force, &temperature);
            let residual_abs = &nondimensional_gibbs_free_energy_per_link - &nondimensional_gibbs_free_energy_per_link_0 - &nondimensional_relative_gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link_0;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
}
mod zero
{
    use super::*;
    use rand::Rng;
    #[test]
    fn relative_gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_gibbs_free_energy_0 = model.relative_gibbs_free_energy(&(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), &temperature);
            assert!(relative_gibbs_free_energy_0.abs() <= BOLTZMANN_CONSTANT*temperature*(number_of_links as f64)*ZERO);
        }
    }
    #[test]
    fn relative_gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_gibbs_free_energy_per_link_0 = model.relative_gibbs_free_energy_per_link(&(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), &temperature);
            assert!(relative_gibbs_free_energy_per_link_0.abs() <= BOLTZMANN_CONSTANT*temperature*ZERO);
        }
    }
    #[test]
    fn nondimensional_relative_gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy_0 = model.nondimensional_relative_gibbs_free_energy(&ZERO, &temperature);
            assert!(nondimensional_relative_gibbs_free_energy_0.abs() <= ZERO);
        }
    }
    #[test]
    fn nondimensional_relative_gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy_per_link_0 = model.nondimensional_relative_gibbs_free_energy_per_link(&ZERO, &temperature);
            assert!(nondimensional_relative_gibbs_free_energy_per_link_0.abs() <= ZERO);
        }
    }
}
mod asymptotic
{
    use super::*;
    use rand::Rng;
    use crate::math::integrate;
    use crate::physics::single_chain::fjc::{
        ZERO,
        POINTS
    };
    #[test]
    fn end_to_end_length()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_link_stiffness|
            {
                let link_stiffness = BOLTZMANN_CONSTANT*temperature/link_length.powi(2)*nondimensional_link_stiffness;
                let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let integrand_numerator = |nondimensional_force: f64|
                {
                    let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
                    (model.end_to_end_length(&force, &temperature) - model.asymptotic.end_to_end_length(&force, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_force: f64|
                {
                    let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
                    model.end_to_end_length(&force, &temperature).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_link_stiffness_medium);
            let residual_rel_2 = residual_rel(parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (-parameters.nondimensional_link_stiffness_medium).exp());
            assert!(residual_rel_2.abs() <= (-parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale).exp());
            assert!(log_log_slope < -2.0);
        }
    }
    #[test]
    fn end_to_end_length_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_link_stiffness|
            {
                let link_stiffness = BOLTZMANN_CONSTANT*temperature/link_length.powi(2)*nondimensional_link_stiffness;
                let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let integrand_numerator = |nondimensional_force: f64|
                {
                    let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
                    (model.end_to_end_length_per_link(&force, &temperature) - model.asymptotic.end_to_end_length_per_link(&force, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_force: f64|
                {
                    let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
                    model.end_to_end_length_per_link(&force, &temperature).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_link_stiffness_medium);
            let residual_rel_2 = residual_rel(parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (-parameters.nondimensional_link_stiffness_medium).exp());
            assert!(residual_rel_2.abs() <= (-parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale).exp());
            assert!(log_log_slope < -2.0);
        }
    }
    #[test]
    fn nondimensional_end_to_end_length()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_link_stiffness|
            {
                let link_stiffness = BOLTZMANN_CONSTANT*temperature/link_length.powi(2)*nondimensional_link_stiffness;
                let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let integrand_numerator = |nondimensional_force: f64|
                {
                    (model.nondimensional_end_to_end_length(&nondimensional_force, &temperature) - model.asymptotic.nondimensional_end_to_end_length(&nondimensional_force, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_force: f64|
                {
                    model.nondimensional_end_to_end_length(&nondimensional_force, &temperature).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_link_stiffness_medium);
            let residual_rel_2 = residual_rel(parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (-parameters.nondimensional_link_stiffness_medium).exp());
            assert!(residual_rel_2.abs() <= (-parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale).exp());
            assert!(log_log_slope < -2.0);
        }
    }
    #[test]
    fn nondimensional_end_to_end_length_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_link_stiffness|
            {
                let link_stiffness = BOLTZMANN_CONSTANT*temperature/link_length.powi(2)*nondimensional_link_stiffness;
                let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let integrand_numerator = |nondimensional_force: f64|
                {
                    (model.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature) - model.asymptotic.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_force: f64|
                {
                    model.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_link_stiffness_medium);
            let residual_rel_2 = residual_rel(parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (-parameters.nondimensional_link_stiffness_medium).exp());
            assert!(residual_rel_2.abs() <= (-parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale).exp());
            assert!(log_log_slope < -2.0);
        }
    }
    #[test]
    fn gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_link_stiffness|
            {
                let link_stiffness = BOLTZMANN_CONSTANT*temperature/link_length.powi(2)*nondimensional_link_stiffness;
                let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let integrand_numerator = |nondimensional_force: f64|
                {
                    let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
                    (model.gibbs_free_energy(&force, &temperature) - model.asymptotic.gibbs_free_energy(&force, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_force: f64|
                {
                    let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
                    model.gibbs_free_energy(&force, &temperature).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_link_stiffness_medium);
            let residual_rel_2 = residual_rel(parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (-parameters.nondimensional_link_stiffness_medium).exp());
            assert!(residual_rel_2.abs() <= (-parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale).exp());
            assert!(log_log_slope < -2.0);
        }
    }
    #[test]
    fn gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_link_stiffness|
            {
                let link_stiffness = BOLTZMANN_CONSTANT*temperature/link_length.powi(2)*nondimensional_link_stiffness;
                let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let integrand_numerator = |nondimensional_force: f64|
                {
                    let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
                    (model.gibbs_free_energy_per_link(&force, &temperature) - model.asymptotic.gibbs_free_energy_per_link(&force, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_force: f64|
                {
                    let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
                    model.gibbs_free_energy_per_link(&force, &temperature).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_link_stiffness_medium);
            let residual_rel_2 = residual_rel(parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (-parameters.nondimensional_link_stiffness_medium).exp());
            assert!(residual_rel_2.abs() <= (-parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale).exp());
            assert!(log_log_slope < -2.0);
        }
    }
    #[test]
    fn relative_gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_link_stiffness|
            {
                let link_stiffness = BOLTZMANN_CONSTANT*temperature/link_length.powi(2)*nondimensional_link_stiffness;
                let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let integrand_numerator = |nondimensional_force: f64|
                {
                    let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
                    (model.relative_gibbs_free_energy(&force, &temperature) - model.asymptotic.relative_gibbs_free_energy(&force, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_force: f64|
                {
                    let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
                    model.relative_gibbs_free_energy(&force, &temperature).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_link_stiffness_medium);
            let residual_rel_2 = residual_rel(parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (-parameters.nondimensional_link_stiffness_medium).exp());
            assert!(residual_rel_2.abs() <= (-parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale).exp());
            assert!(log_log_slope < -2.0);
        }
    }
    #[test]
    fn relative_gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_link_stiffness|
            {
                let link_stiffness = BOLTZMANN_CONSTANT*temperature/link_length.powi(2)*nondimensional_link_stiffness;
                let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let integrand_numerator = |nondimensional_force: f64|
                {
                    let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
                    (model.relative_gibbs_free_energy_per_link(&force, &temperature) - model.asymptotic.relative_gibbs_free_energy_per_link(&force, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_force: f64|
                {
                    let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
                    model.relative_gibbs_free_energy_per_link(&force, &temperature).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_link_stiffness_medium);
            let residual_rel_2 = residual_rel(parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (-parameters.nondimensional_link_stiffness_medium).exp());
            assert!(residual_rel_2.abs() <= (-parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale).exp());
            assert!(log_log_slope < -2.0);
        }
    }
    #[test]
    fn nondimensional_gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_link_stiffness|
            {
                let link_stiffness = BOLTZMANN_CONSTANT*temperature/link_length.powi(2)*nondimensional_link_stiffness;
                let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let integrand_numerator = |nondimensional_force: f64|
                {
                    (model.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature) - model.asymptotic.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_force: f64|
                {
                    model.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_link_stiffness_medium);
            let residual_rel_2 = residual_rel(parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (-parameters.nondimensional_link_stiffness_medium).exp());
            assert!(residual_rel_2.abs() <= (-parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale).exp());
            assert!(log_log_slope < -2.0);
        }
    }
    #[test]
    fn nondimensional_gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_link_stiffness|
            {
                let link_stiffness = BOLTZMANN_CONSTANT*temperature/link_length.powi(2)*nondimensional_link_stiffness;
                let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let integrand_numerator = |nondimensional_force: f64|
                {
                    (model.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature) - model.asymptotic.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_force: f64|
                {
                    model.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_link_stiffness_medium);
            let residual_rel_2 = residual_rel(parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (-parameters.nondimensional_link_stiffness_medium).exp());
            assert!(residual_rel_2.abs() <= (-parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale).exp());
            assert!(log_log_slope < -2.0);
        }
    }
    #[test]
    fn nondimensional_relative_gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_link_stiffness|
            {
                let link_stiffness = BOLTZMANN_CONSTANT*temperature/link_length.powi(2)*nondimensional_link_stiffness;
                let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let integrand_numerator = |nondimensional_force: f64|
                {
                    (model.nondimensional_relative_gibbs_free_energy(&nondimensional_force, &temperature) - model.asymptotic.nondimensional_relative_gibbs_free_energy(&nondimensional_force, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_force: f64|
                {
                    model.nondimensional_relative_gibbs_free_energy(&nondimensional_force, &temperature).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_link_stiffness_medium);
            let residual_rel_2 = residual_rel(parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (-parameters.nondimensional_link_stiffness_medium).exp());
            assert!(residual_rel_2.abs() <= (-parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale).exp());
            assert!(log_log_slope < -2.0);
        }
    }
    #[test]
    fn nondimensional_relative_gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_link_stiffness|
            {
                let link_stiffness = BOLTZMANN_CONSTANT*temperature/link_length.powi(2)*nondimensional_link_stiffness;
                let model = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let integrand_numerator = |nondimensional_force: f64|
                {
                    (model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force, &temperature) - model.asymptotic.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_force: f64|
                {
                    model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force, &temperature).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_force_scale, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_link_stiffness_medium);
            let residual_rel_2 = residual_rel(parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (-parameters.nondimensional_link_stiffness_medium).exp());
            assert!(residual_rel_2.abs() <= (-parameters.nondimensional_link_stiffness_medium*parameters.log_log_scale).exp());
            assert!(log_log_slope < -2.0);
        }
    }
}
mod implementations
{
    use super::*;
    mod asymptotic
    {
        use super::*;
        #[test]
        fn access()
        {
            let parameters = Parameters::default();
            let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).asymptotic.nondimensional_end_to_end_length_per_link(&parameters.nondimensional_force_reference, &parameters.temperature_reference);
        }
        mod legendre
        {
            use super::*;
            #[test]
            fn access()
            {
                let parameters = Parameters::default();
                let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).asymptotic.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&parameters.nondimensional_force_reference, &parameters.temperature_reference);
            }
        }
        mod alternative
        {
            use super::*;
            #[test]
            fn access()
            {
                let parameters = Parameters::default();
                let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).asymptotic.alternative.nondimensional_end_to_end_length_per_link(&parameters.nondimensional_force_reference, &parameters.temperature_reference);
            }
            mod legendre
            {
                use super::*;
                #[test]
                fn access()
                {
                    let parameters = Parameters::default();
                    let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).asymptotic.alternative.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&parameters.nondimensional_force_reference, &parameters.temperature_reference);
                }
            }
        }
        mod reduced
        {
            use super::*;
            #[test]
            fn access()
            {
                let parameters = Parameters::default();
                let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).asymptotic.reduced.nondimensional_end_to_end_length_per_link(&parameters.nondimensional_force_reference, &parameters.temperature_reference);
            }
            mod legendre
            {
                use super::*;
                #[test]
                fn access()
                {
                    let parameters = Parameters::default();
                    let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).asymptotic.reduced.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&parameters.nondimensional_force_reference, &parameters.temperature_reference);
                }
            }
        }
    }
}