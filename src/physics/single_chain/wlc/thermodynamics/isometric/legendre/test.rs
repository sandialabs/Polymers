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
        let _ = WLC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.persistance_length_reference);
    }
    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            assert_eq!(number_of_links, WLC::init(number_of_links, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.persistance_length_reference).number_of_links);
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
            assert_eq!(link_length, WLC::init(parameters.number_of_links_minimum, link_length, parameters.hinge_mass_reference, parameters.persistance_length_reference).link_length);
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
            assert_eq!(hinge_mass, WLC::init(parameters.number_of_links_minimum, parameters.link_length_reference, hinge_mass, parameters.persistance_length_reference).hinge_mass);
        }
    }
    #[test]
    fn persistance_length()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(persistance_length, WLC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, persistance_length).persistance_length);
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
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            assert_eq!(number_of_links, model.number_of_links);
            assert_eq!(link_length, model.link_length);
            assert_eq!(hinge_mass, model.hinge_mass);
            assert_eq!(persistance_length, model.persistance_length);
        }
    }
}
mod nondimensional
{
    use super::*;
    use rand::Rng;
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let gibbs_free_energy = model.gibbs_free_energy(&end_to_end_length, &temperature);
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&end_to_end_length, &temperature);
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_end_to_end_length_per_link);
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&end_to_end_length, &temperature);
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&end_to_end_length, &temperature);
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
    fn gibbs_free_energy()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let gibbs_free_energy = model.gibbs_free_energy(&end_to_end_length, &temperature);
            let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&end_to_end_length, &temperature);
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&end_to_end_length, &temperature);
            let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&end_to_end_length, &temperature);
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_end_to_end_length_per_link);
            let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
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
    use crate::physics::single_chain::ZERO;
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let gibbs_free_energy = model.gibbs_free_energy(&end_to_end_length, &temperature);
            let gibbs_free_energy_0 = model.gibbs_free_energy(&(ZERO*(number_of_links as f64)*link_length), &temperature);
            let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&end_to_end_length, &temperature);
            let residual_abs = &gibbs_free_energy - &gibbs_free_energy_0 - &relative_gibbs_free_energy;
            let residual_rel = &residual_abs/&gibbs_free_energy_0;
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&end_to_end_length, &temperature);
            let gibbs_free_energy_per_link_0 = model.gibbs_free_energy_per_link(&(ZERO*(number_of_links as f64)*link_length), &temperature);
            let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&end_to_end_length, &temperature);
            let residual_abs = &gibbs_free_energy_per_link - &gibbs_free_energy_per_link_0 - &relative_gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&gibbs_free_energy_per_link_0;
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let nondimensional_gibbs_free_energy_0 = model.nondimensional_gibbs_free_energy(&ZERO, &temperature);
            let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_end_to_end_length_per_link);
            let residual_abs = &nondimensional_gibbs_free_energy - &nondimensional_gibbs_free_energy_0 - &nondimensional_relative_gibbs_free_energy;
            let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_0;
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
            let nondimensional_gibbs_free_energy_per_link_0 = model.nondimensional_gibbs_free_energy_per_link(&ZERO, &temperature);
            let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
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
    use crate::physics::single_chain::ZERO;
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_gibbs_free_energy_0 = model.relative_gibbs_free_energy(&(ZERO*(number_of_links as f64)*link_length), &temperature);
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_gibbs_free_energy_per_link_0 = model.relative_gibbs_free_energy_per_link(&(ZERO*(number_of_links as f64)*link_length), &temperature);
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_relative_gibbs_free_energy_0 = model.nondimensional_relative_gibbs_free_energy(&ZERO);
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
            let persistance_length = parameters.persistance_length_reference + parameters.persistance_length_scale*(0.5 - rng.gen::<f64>());
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_relative_gibbs_free_energy_per_link_0 = model.nondimensional_relative_gibbs_free_energy_per_link(&ZERO);
            assert!(nondimensional_relative_gibbs_free_energy_per_link_0.abs() <= ZERO);
        }
    }
}