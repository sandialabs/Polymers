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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let end_to_end_length = model.end_to_end_length(&potential_distance, &potential_stiffness, &temperature);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let end_to_end_length = model.end_to_end_length_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = &end_to_end_length/link_length - &nondimensional_end_to_end_length;
            let residual_rel = &residual_abs/&nondimensional_end_to_end_length;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn force()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force = model.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let force = model.force(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = &force/BOLTZMANN_CONSTANT/temperature*link_length - &nondimensional_force;
            let residual_rel = &residual_abs/&nondimensional_force;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let helmholtz_free_energy = model.helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = helmholtz_free_energy/BOLTZMANN_CONSTANT/temperature - nondimensional_helmholtz_free_energy;
            let residual_rel = residual_abs/nondimensional_helmholtz_free_energy;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = helmholtz_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - nondimensional_helmholtz_free_energy_per_link;
            let residual_rel = residual_abs/nondimensional_helmholtz_free_energy_per_link;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = relative_helmholtz_free_energy/BOLTZMANN_CONSTANT/temperature - nondimensional_relative_helmholtz_free_energy;
            let residual_rel = residual_abs/nondimensional_relative_helmholtz_free_energy;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = relative_helmholtz_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - nondimensional_relative_helmholtz_free_energy_per_link;
            let residual_rel = residual_abs/nondimensional_relative_helmholtz_free_energy_per_link;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let gibbs_free_energy = model.gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = gibbs_free_energy/BOLTZMANN_CONSTANT/temperature - nondimensional_gibbs_free_energy;
            let residual_rel = residual_abs/nondimensional_gibbs_free_energy;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = gibbs_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - nondimensional_gibbs_free_energy_per_link;
            let residual_rel = residual_abs/nondimensional_gibbs_free_energy_per_link;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = relative_gibbs_free_energy/BOLTZMANN_CONSTANT/temperature - nondimensional_relative_gibbs_free_energy;
            let residual_rel = residual_abs/nondimensional_relative_gibbs_free_energy;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = relative_gibbs_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - nondimensional_relative_gibbs_free_energy_per_link;
            let residual_rel = residual_abs/nondimensional_relative_gibbs_free_energy_per_link;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let end_to_end_length = model.end_to_end_length(&potential_distance, &potential_stiffness, &temperature);
            let end_to_end_length_per_link = model.end_to_end_length_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = end_to_end_length/(number_of_links as f64) - end_to_end_length_per_link;
            let residual_rel = residual_abs/end_to_end_length_per_link;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let residual_abs = nondimensional_end_to_end_length/(number_of_links as f64) - nondimensional_end_to_end_length_per_link;
            let residual_rel = residual_abs/nondimensional_end_to_end_length_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let helmholtz_free_energy = model.helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = helmholtz_free_energy/(number_of_links as f64) - helmholtz_free_energy_per_link;
            let residual_rel = residual_abs/helmholtz_free_energy_per_link;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = relative_helmholtz_free_energy/(number_of_links as f64) - relative_helmholtz_free_energy_per_link;
            let residual_rel = residual_abs/relative_helmholtz_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn nondimensional_helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
            let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
            let residual_abs = nondimensional_helmholtz_free_energy/(number_of_links as f64) - nondimensional_helmholtz_free_energy_per_link;
            let residual_rel = residual_abs/nondimensional_helmholtz_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let residual_abs = nondimensional_relative_helmholtz_free_energy/(number_of_links as f64) - nondimensional_relative_helmholtz_free_energy_per_link;
            let residual_rel = residual_abs/nondimensional_relative_helmholtz_free_energy_per_link;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let gibbs_free_energy = model.gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = gibbs_free_energy/(number_of_links as f64) - gibbs_free_energy_per_link;
            let residual_rel = residual_abs/gibbs_free_energy_per_link;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = relative_gibbs_free_energy/(number_of_links as f64) - relative_gibbs_free_energy_per_link;
            let residual_rel = residual_abs/relative_gibbs_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
            let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
            let residual_abs = nondimensional_gibbs_free_energy/(number_of_links as f64) - nondimensional_gibbs_free_energy_per_link;
            let residual_rel = residual_abs/nondimensional_gibbs_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let residual_abs = nondimensional_relative_gibbs_free_energy/(number_of_links as f64) - nondimensional_relative_gibbs_free_energy_per_link;
            let residual_rel = residual_abs/nondimensional_relative_gibbs_free_energy_per_link;
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
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let helmholtz_free_energy = model.helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let helmholtz_free_energy_0 = model.helmholtz_free_energy(&ZERO, &potential_stiffness, &temperature);
            let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = &helmholtz_free_energy - &helmholtz_free_energy_0 - &relative_helmholtz_free_energy;
            let residual_rel = &residual_abs/&helmholtz_free_energy_0;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let helmholtz_free_energy_per_link_0 = model.helmholtz_free_energy_per_link(&ZERO, &potential_stiffness, &temperature);
            let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = &helmholtz_free_energy_per_link - &helmholtz_free_energy_per_link_0 - &relative_helmholtz_free_energy_per_link;
            let residual_rel = &residual_abs/&helmholtz_free_energy_per_link_0;
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn nondimensional_helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
            let nondimensional_helmholtz_free_energy_0 = model.nondimensional_helmholtz_free_energy(&ZERO, &nondimensional_potential_stiffness, &temperature);
            let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let residual_abs = &nondimensional_helmholtz_free_energy - &nondimensional_helmholtz_free_energy_0 - &nondimensional_relative_helmholtz_free_energy;
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_0;
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn nondimensional_helmholtz_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
            let nondimensional_helmholtz_free_energy_per_link_0 = model.nondimensional_helmholtz_free_energy_per_link(&ZERO, &nondimensional_potential_stiffness, &temperature);
            let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let residual_abs = &nondimensional_helmholtz_free_energy_per_link - &nondimensional_helmholtz_free_energy_per_link_0 - &nondimensional_relative_helmholtz_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link_0;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let gibbs_free_energy = model.gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let gibbs_free_energy_0 = model.gibbs_free_energy(&ZERO, &potential_stiffness, &temperature);
            let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = &gibbs_free_energy - &gibbs_free_energy_0 - &relative_gibbs_free_energy;
            let residual_rel = &residual_abs/&gibbs_free_energy_0;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let gibbs_free_energy_per_link_0 = model.gibbs_free_energy_per_link(&ZERO, &potential_stiffness, &temperature);
            let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = &gibbs_free_energy_per_link - &gibbs_free_energy_per_link_0 - &relative_gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&gibbs_free_energy_per_link_0;
            assert!(residual_rel.abs() <= parameters.rel_tol);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
            let nondimensional_gibbs_free_energy_0 = model.nondimensional_gibbs_free_energy(&ZERO, &nondimensional_potential_stiffness, &temperature);
            let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let residual_abs = &nondimensional_gibbs_free_energy - &nondimensional_gibbs_free_energy_0 - &nondimensional_relative_gibbs_free_energy;
            let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_0;
            assert!(residual_rel.abs() <= parameters.rel_tol);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
            let nondimensional_gibbs_free_energy_per_link_0 = model.nondimensional_gibbs_free_energy_per_link(&ZERO, &nondimensional_potential_stiffness, &temperature);
            let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let residual_abs = &nondimensional_gibbs_free_energy_per_link - &nondimensional_gibbs_free_energy_per_link_0 - &nondimensional_relative_gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link_0;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let relative_helmholtz_free_energy_0 = model.relative_helmholtz_free_energy(&(ZERO*(number_of_links as f64)*link_length), &potential_stiffness, &temperature);
            assert!(relative_helmholtz_free_energy_0.abs() <= BOLTZMANN_CONSTANT*temperature*(number_of_links as f64)*ZERO);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let relative_helmholtz_free_energy_per_link_0 = model.relative_helmholtz_free_energy_per_link(&(ZERO*(number_of_links as f64)*link_length), &potential_stiffness, &temperature);
            assert!(relative_helmholtz_free_energy_per_link_0.abs() <= BOLTZMANN_CONSTANT*temperature*ZERO);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_helmholtz_free_energy_0 = model.nondimensional_relative_helmholtz_free_energy(&ZERO, &nondimensional_potential_stiffness);
            assert!(nondimensional_relative_helmholtz_free_energy_0.abs() <= (number_of_links as f64)*ZERO);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_helmholtz_free_energy_per_link_0 = model.nondimensional_relative_helmholtz_free_energy_per_link(&ZERO, &nondimensional_potential_stiffness);
            assert!(nondimensional_relative_helmholtz_free_energy_per_link_0.abs() <= ZERO);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let relative_gibbs_free_energy_0 = model.relative_gibbs_free_energy(&ZERO, &potential_stiffness, &temperature);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let relative_gibbs_free_energy_per_link_0 = model.relative_gibbs_free_energy_per_link(&ZERO, &potential_stiffness, &temperature);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy_0 = model.nondimensional_relative_gibbs_free_energy(&ZERO, &nondimensional_potential_stiffness);
            assert!(nondimensional_relative_gibbs_free_energy_0.abs() <= (number_of_links as f64)*ZERO);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy_per_link_0 = model.nondimensional_relative_gibbs_free_energy_per_link(&ZERO, &nondimensional_potential_stiffness);
            assert!(nondimensional_relative_gibbs_free_energy_per_link_0.abs() <= ZERO);
        }
    }
}
mod strong_potential
{
    use super::*;
    use rand::Rng;
    use crate::physics::single_chain::
    {
        ZERO,
        POINTS,
        test::integrate
    };
    #[test]
    fn force()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_numerator = |potential_distance: f64|
                {
                    (model.force(&potential_distance, &potential_stiffness, &temperature) - model.asymptotic.strong_potential.force(&potential_distance, &potential_stiffness, &temperature)).powi(2)
                };
                let integrand_denominator = |potential_distance: f64|
                {
                    let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                    (model.force(&potential_distance, &potential_stiffness, &temperature)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                let denominator = integrate(integrand_denominator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn nondimensional_force()
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
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.asymptotic.strong_potential.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_numerator = |potential_distance: f64|
                {
                    (model.helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature) - model.helmholtz_free_energy(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature) - model.asymptotic.strong_potential.helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature) + model.asymptotic.strong_potential.helmholtz_free_energy(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature)).powi(2)
                };
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_denominator = |potential_distance: f64|
                {
                    (model.helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature) - model.helmholtz_free_energy(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                let denominator = integrate(integrand_denominator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn helmholtz_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_numerator = |potential_distance: f64|
                {
                    (model.helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) - model.helmholtz_free_energy_per_link(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature) - model.asymptotic.strong_potential.helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) + model.asymptotic.strong_potential.helmholtz_free_energy_per_link(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature)).powi(2)
                };
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_denominator = |potential_distance: f64|
                {
                    (model.helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) - model.helmholtz_free_energy_per_link(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                let denominator = integrate(integrand_denominator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_numerator = |potential_distance: f64|
                {
                    (model.relative_helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature) - model.relative_helmholtz_free_energy(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature) - model.asymptotic.strong_potential.relative_helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature) + model.asymptotic.strong_potential.relative_helmholtz_free_energy(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature)).powi(2)
                };
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_denominator = |potential_distance: f64|
                {
                    (model.relative_helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature) - model.relative_helmholtz_free_energy(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                let denominator = integrate(integrand_denominator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_numerator = |potential_distance: f64|
                {
                    (model.relative_helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) - model.relative_helmholtz_free_energy_per_link(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature) - model.asymptotic.strong_potential.relative_helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) + model.asymptotic.strong_potential.relative_helmholtz_free_energy_per_link(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature)).powi(2)
                };
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_denominator = |potential_distance: f64|
                {
                    (model.relative_helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) - model.relative_helmholtz_free_energy_per_link(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                let denominator = integrate(integrand_denominator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn nondimensional_helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) - model.nondimensional_helmholtz_free_energy(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness, &temperature) - model.asymptotic.strong_potential.nondimensional_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) + model.asymptotic.strong_potential.nondimensional_helmholtz_free_energy(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) - model.nondimensional_helmholtz_free_energy(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness, &temperature)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn nondimensional_helmholtz_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) - model.nondimensional_helmholtz_free_energy_per_link(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness, &temperature) - model.asymptotic.strong_potential.nondimensional_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) + model.asymptotic.strong_potential.nondimensional_helmholtz_free_energy_per_link(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) - model.nondimensional_helmholtz_free_energy_per_link(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness, &temperature)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
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
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_relative_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.nondimensional_relative_helmholtz_free_energy(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness) - model.asymptotic.strong_potential.nondimensional_relative_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness) + model.asymptotic.strong_potential.nondimensional_relative_helmholtz_free_energy(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_relative_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.nondimensional_relative_helmholtz_free_energy(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
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
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.nondimensional_relative_helmholtz_free_energy_per_link(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness) - model.asymptotic.strong_potential.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness) + model.asymptotic.strong_potential.nondimensional_relative_helmholtz_free_energy_per_link(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.nondimensional_relative_helmholtz_free_energy_per_link(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
            let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
}
mod weak_potential
{
    use super::*;
    use rand::Rng;
    use crate::physics::single_chain::test::integrate;
    #[test]
    fn end_to_end_length()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                    (model.end_to_end_length(&potential_distance, &potential_stiffness, &temperature) - model.asymptotic.weak_potential.end_to_end_length(&potential_distance, &potential_stiffness, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                    model.end_to_end_length(&potential_distance, &potential_stiffness, &temperature).powi(2)
                };
                let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
            // let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (parameters.nondimensional_potential_stiffness_small).powi(2));
            assert!(residual_rel_2.abs() <= (parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale).powi(2));
            // assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn end_to_end_length_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                    (model.end_to_end_length_per_link(&potential_distance, &potential_stiffness, &temperature) - model.asymptotic.weak_potential.end_to_end_length_per_link(&potential_distance, &potential_stiffness, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                    model.end_to_end_length_per_link(&potential_distance, &potential_stiffness, &temperature).powi(2)
                };
                let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
            // let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (parameters.nondimensional_potential_stiffness_small).powi(2));
            assert!(residual_rel_2.abs() <= (parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale).powi(2));
            // assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn nondimensional_end_to_end_length()
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
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_end_to_end_length(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.asymptotic.weak_potential.nondimensional_end_to_end_length(&nondimensional_potential_distance, &nondimensional_potential_stiffness)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    model.nondimensional_end_to_end_length(&nondimensional_potential_distance, &nondimensional_potential_stiffness).powi(2)
                };
                let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
            // let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (parameters.nondimensional_potential_stiffness_small).powi(2));
            assert!(residual_rel_2.abs() <= (parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale).powi(2));
            // assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn nondimensional_end_to_end_length_per_link()
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
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_end_to_end_length_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.asymptotic.weak_potential.nondimensional_end_to_end_length_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    model.nondimensional_end_to_end_length_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness).powi(2)
                };
                let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
            // let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (parameters.nondimensional_potential_stiffness_small).powi(2));
            assert!(residual_rel_2.abs() <= (parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale).powi(2));
            // assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                    (model.gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature) - model.gibbs_free_energy(&((number_of_links as f64)*link_length*parameters.nondimensional_potential_distance_large_1), &potential_stiffness, &temperature) - model.asymptotic.weak_potential.gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature) + model.asymptotic.weak_potential.gibbs_free_energy(&((number_of_links as f64)*link_length*parameters.nondimensional_potential_distance_large_1), &potential_stiffness, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                    (model.gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature) - model.gibbs_free_energy(&((number_of_links as f64)*link_length*parameters.nondimensional_potential_distance_large_1), &potential_stiffness, &temperature)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
            // let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (parameters.nondimensional_potential_stiffness_small).powi(2));
            assert!(residual_rel_2.abs() <= (parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale).powi(2));
            // assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                    (model.gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) - model.gibbs_free_energy_per_link(&((number_of_links as f64)*link_length*parameters.nondimensional_potential_distance_large_1), &potential_stiffness, &temperature) - model.asymptotic.weak_potential.gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) + model.asymptotic.weak_potential.gibbs_free_energy_per_link(&((number_of_links as f64)*link_length*parameters.nondimensional_potential_distance_large_1), &potential_stiffness, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                    (model.gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) - model.gibbs_free_energy_per_link(&((number_of_links as f64)*link_length*parameters.nondimensional_potential_distance_large_1), &potential_stiffness, &temperature)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
            // let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (parameters.nondimensional_potential_stiffness_small).powi(2));
            assert!(residual_rel_2.abs() <= (parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale).powi(2));
            // assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn relative_gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                    (model.relative_gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature) - model.relative_gibbs_free_energy(&((number_of_links as f64)*link_length*parameters.nondimensional_potential_distance_large_1), &potential_stiffness, &temperature) - model.asymptotic.weak_potential.relative_gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature) + model.asymptotic.weak_potential.relative_gibbs_free_energy(&((number_of_links as f64)*link_length*parameters.nondimensional_potential_distance_large_1), &potential_stiffness, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                    (model.relative_gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature) - model.relative_gibbs_free_energy(&((number_of_links as f64)*link_length*parameters.nondimensional_potential_distance_large_1), &potential_stiffness, &temperature)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
            // let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (parameters.nondimensional_potential_stiffness_small).powi(2));
            assert!(residual_rel_2.abs() <= (parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale).powi(2));
            // assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn relative_gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                    (model.relative_gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) - model.relative_gibbs_free_energy_per_link(&((number_of_links as f64)*link_length*parameters.nondimensional_potential_distance_large_1), &potential_stiffness, &temperature) - model.asymptotic.weak_potential.relative_gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) + model.asymptotic.weak_potential.relative_gibbs_free_energy_per_link(&((number_of_links as f64)*link_length*parameters.nondimensional_potential_distance_large_1), &potential_stiffness, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                    (model.relative_gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) - model.relative_gibbs_free_energy_per_link(&((number_of_links as f64)*link_length*parameters.nondimensional_potential_distance_large_1), &potential_stiffness, &temperature)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
            // let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (parameters.nondimensional_potential_stiffness_small).powi(2));
            assert!(residual_rel_2.abs() <= (parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale).powi(2));
            // assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn nondimensional_gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) - model.nondimensional_gibbs_free_energy(&parameters.nondimensional_potential_distance_large_1, &nondimensional_potential_stiffness, &temperature) - model.asymptotic.weak_potential.nondimensional_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) + model.asymptotic.weak_potential.nondimensional_gibbs_free_energy(&parameters.nondimensional_potential_distance_large_1, &nondimensional_potential_stiffness, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) - model.nondimensional_gibbs_free_energy(&parameters.nondimensional_potential_distance_large_1, &nondimensional_potential_stiffness, &temperature)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
            // let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (parameters.nondimensional_potential_stiffness_small).powi(2));
            assert!(residual_rel_2.abs() <= (parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale).powi(2));
            // assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn nondimensional_gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum - parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let residual_rel = |nondimensional_potential_stiffness|
            {
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) - model.nondimensional_gibbs_free_energy_per_link(&parameters.nondimensional_potential_distance_large_1, &nondimensional_potential_stiffness, &temperature) - model.asymptotic.weak_potential.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) + model.asymptotic.weak_potential.nondimensional_gibbs_free_energy_per_link(&parameters.nondimensional_potential_distance_large_1, &nondimensional_potential_stiffness, &temperature)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) - model.nondimensional_gibbs_free_energy_per_link(&parameters.nondimensional_potential_distance_large_1, &nondimensional_potential_stiffness, &temperature)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
            // let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (parameters.nondimensional_potential_stiffness_small).powi(2));
            assert!(residual_rel_2.abs() <= (parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale).powi(2));
            // assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn nondimensional_relative_gibbs_free_energy()
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
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.nondimensional_relative_gibbs_free_energy(&parameters.nondimensional_potential_distance_large_1, &nondimensional_potential_stiffness) - model.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness) + model.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy(&parameters.nondimensional_potential_distance_large_1, &nondimensional_potential_stiffness)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.nondimensional_relative_gibbs_free_energy(&parameters.nondimensional_potential_distance_large_1, &nondimensional_potential_stiffness)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
            // let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (parameters.nondimensional_potential_stiffness_small).powi(2));
            assert!(residual_rel_2.abs() <= (parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale).powi(2));
            // assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
    #[test]
    fn nondimensional_relative_gibbs_free_energy_per_link()
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
                let integrand_numerator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.nondimensional_relative_gibbs_free_energy_per_link(&parameters.nondimensional_potential_distance_large_1, &nondimensional_potential_stiffness) - model.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness) + model.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy_per_link(&parameters.nondimensional_potential_distance_large_1, &nondimensional_potential_stiffness)).powi(2)
                };
                let integrand_denominator = |nondimensional_potential_distance: f64|
                {
                    (model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.nondimensional_relative_gibbs_free_energy_per_link(&parameters.nondimensional_potential_distance_large_1, &nondimensional_potential_stiffness)).powi(2)
                };
                let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                (numerator/denominator).sqrt()
            };
            let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
            let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
            // let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
            assert!(residual_rel_1.abs() <= (parameters.nondimensional_potential_stiffness_small).powi(2));
            assert!(residual_rel_2.abs() <= (parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale).powi(2));
            // assert!((log_log_slope/2.0 + 1.0).abs() <= parameters.log_log_tol);
        }
    }
}
