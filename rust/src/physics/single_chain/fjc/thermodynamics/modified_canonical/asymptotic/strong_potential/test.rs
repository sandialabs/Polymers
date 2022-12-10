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
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_small*(0.5 + (0.5 - rng.gen::<f64>()));
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force = model.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
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
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_small*(0.5 + (0.5 - rng.gen::<f64>()));
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
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
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_small*(0.5 + (0.5 - rng.gen::<f64>()));
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = helmholtz_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - nondimensional_helmholtz_free_energy_per_link;
            let residual_rel = residual_abs/nondimensional_helmholtz_free_energy_per_link;
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
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_small*(0.5 + (0.5 - rng.gen::<f64>()));
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
            let helmholtz_free_energy = model.helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = helmholtz_free_energy/(number_of_links as f64) - helmholtz_free_energy_per_link;
            let residual_rel = residual_abs/helmholtz_free_energy_per_link;
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
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_small*(0.5 + (0.5 - rng.gen::<f64>()));
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
}