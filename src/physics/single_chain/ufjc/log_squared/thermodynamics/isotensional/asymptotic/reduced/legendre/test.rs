#![cfg(test)]
use super::*;
use crate::physics::single_chain::test::Parameters;
mod base
{
    use super::*;
    use rand::Rng;
    #[test]
    fn init()
    {
        let parameters = Parameters::default();
        let _ = LOGSQUAREDFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference);
    }
    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            assert_eq!(number_of_links, LOGSQUAREDFJC::init(number_of_links, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).number_of_links);
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
            assert_eq!(link_length, LOGSQUAREDFJC::init(parameters.number_of_links_minimum, link_length, parameters.hinge_mass_reference, parameters.link_stiffness_reference).link_length);
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
            assert_eq!(hinge_mass, LOGSQUAREDFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, hinge_mass, parameters.link_stiffness_reference).hinge_mass);
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
            assert_eq!(link_stiffness, LOGSQUAREDFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, link_stiffness).link_stiffness);
        }
    }
    #[test]
    fn all_parameters()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_length = rng.gen::<f64>();
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            assert_eq!(number_of_links, model.number_of_links);
            assert_eq!(link_length, model.link_length);
            assert_eq!(hinge_mass, model.hinge_mass);
            assert_eq!(link_stiffness, model.link_stiffness);
        }
    }
}
mod nondimensional
{
    use super::*;
    use rand::Rng;
    #[test]
    fn helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/1.0_f64.exp();
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_force, &temperature);
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let helmholtz_free_energy = model.helmholtz_free_energy(&force, &temperature);
            let residual_abs = &helmholtz_free_energy/BOLTZMANN_CONSTANT/temperature - &nondimensional_helmholtz_free_energy;
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn helmholtz_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/1.0_f64.exp();
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_force, &temperature);
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&force, &temperature);
            let residual_abs = &helmholtz_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_helmholtz_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn relative_helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/1.0_f64.exp();
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_force, &temperature);
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&force, &temperature);
            let residual_abs = &relative_helmholtz_free_energy/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_helmholtz_free_energy;
            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn relative_helmholtz_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/1.0_f64.exp();
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force, &temperature);
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&force, &temperature);
            let residual_abs = &relative_helmholtz_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_helmholtz_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
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
    fn helmholtz_free_energy()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/1.0_f64.exp();
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let helmholtz_free_energy = model.helmholtz_free_energy(&force, &temperature);
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&force, &temperature);
            let residual_abs = &helmholtz_free_energy/(number_of_links as f64) - &helmholtz_free_energy_per_link;
            let residual_rel = &residual_abs/&helmholtz_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn relative_helmholtz_free_energy()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/1.0_f64.exp();
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&force, &temperature);
            let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&force, &temperature);
            let residual_abs = &relative_helmholtz_free_energy/(number_of_links as f64) - &relative_helmholtz_free_energy_per_link;
            let residual_rel = &residual_abs/&relative_helmholtz_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn nondimensional_helmholtz_free_energy()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/1.0_f64.exp();
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_force, &temperature);
            let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_force, &temperature);
            let residual_abs = &nondimensional_helmholtz_free_energy/(number_of_links as f64) - &nondimensional_helmholtz_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn nondimensional_relative_helmholtz_free_energy()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/1.0_f64.exp();
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_force, &temperature);
            let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force, &temperature);
            let residual_abs = &nondimensional_relative_helmholtz_free_energy/(number_of_links as f64) - &nondimensional_relative_helmholtz_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
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
    fn helmholtz_free_energy()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/1.0_f64.exp();
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let helmholtz_free_energy = model.helmholtz_free_energy(&force, &temperature);
            let helmholtz_free_energy_0 = model.helmholtz_free_energy(&(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), &temperature);
            let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&force, &temperature);
            let residual_abs = &helmholtz_free_energy - &helmholtz_free_energy_0 - &relative_helmholtz_free_energy;
            let residual_rel = &residual_abs/&helmholtz_free_energy_0;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn helmholtz_free_energy_per_link()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/1.0_f64.exp();
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&force, &temperature);
            let helmholtz_free_energy_per_link_0 = model.helmholtz_free_energy_per_link(&(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), &temperature);
            let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&force, &temperature);
            let residual_abs = &helmholtz_free_energy_per_link - &helmholtz_free_energy_per_link_0 - &relative_helmholtz_free_energy_per_link;
            let residual_rel = &residual_abs/&helmholtz_free_energy_per_link_0;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn nondimensional_helmholtz_free_energy()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/1.0_f64.exp();
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_force, &temperature);
            let nondimensional_helmholtz_free_energy_0 = model.nondimensional_helmholtz_free_energy(&ZERO, &temperature);
            let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_force, &temperature);
            let residual_abs = &nondimensional_helmholtz_free_energy - &nondimensional_helmholtz_free_energy_0 - &nondimensional_relative_helmholtz_free_energy;
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_0;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn nondimensional_helmholtz_free_energy_per_link()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/1.0_f64.exp();
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_force, &temperature);
            let nondimensional_helmholtz_free_energy_per_link_0 = model.nondimensional_helmholtz_free_energy_per_link(&ZERO, &temperature);
            let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force, &temperature);
            let residual_abs = &nondimensional_helmholtz_free_energy_per_link - &nondimensional_helmholtz_free_energy_per_link_0 - &nondimensional_relative_helmholtz_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link_0;
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
    fn relative_helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_helmholtz_free_energy_0 = model.relative_helmholtz_free_energy(&(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), &temperature);
            assert!(relative_helmholtz_free_energy_0.abs() <= ZERO);
        }
    }
    #[test]
    fn relative_helmholtz_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_helmholtz_free_energy_per_link_0 = model.relative_helmholtz_free_energy_per_link(&(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), &temperature);
            assert!(relative_helmholtz_free_energy_per_link_0.abs() <= ZERO);
        }
    }
    #[test]
    fn nondimensional_relative_helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_helmholtz_free_energy_0 = model.nondimensional_relative_helmholtz_free_energy(&ZERO, &temperature);
            assert!(nondimensional_relative_helmholtz_free_energy_0.abs() <= ZERO);
        }
    }
    #[test]
    fn nondimensional_relative_helmholtz_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LOGSQUAREDFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_helmholtz_free_energy_per_link_0 = model.nondimensional_relative_helmholtz_free_energy_per_link(&ZERO, &temperature);
            assert!(nondimensional_relative_helmholtz_free_energy_per_link_0.abs() <= ZERO);
        }
    }
}
