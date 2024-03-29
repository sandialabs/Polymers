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
        let _ = Ideal::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference);
    }
    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            assert_eq!(number_of_links, Ideal::init(number_of_links, parameters.link_length_reference, parameters.hinge_mass_reference).number_of_links);
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
            assert_eq!(link_length, Ideal::init(parameters.number_of_links_minimum, link_length, parameters.hinge_mass_reference).link_length);
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
            assert_eq!(hinge_mass, Ideal::init(parameters.number_of_links_minimum, parameters.link_length_reference, hinge_mass).hinge_mass);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            assert_eq!(number_of_links, model.number_of_links);
            assert_eq!(link_length, model.link_length);
            assert_eq!(hinge_mass, model.hinge_mass);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_force);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&nondimensional_force);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_force);
            let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&nondimensional_force);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
            let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let gibbs_free_energy = model.gibbs_free_energy(&force, &temperature);
            let gibbs_free_energy_0 = model.gibbs_free_energy(&0.0, &temperature);
            let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&force, &temperature);
            let residual_abs = &gibbs_free_energy - &gibbs_free_energy_0 - &relative_gibbs_free_energy;
            let residual_rel = &residual_abs/&gibbs_free_energy;
            assert!(residual_abs.abs() <= BOLTZMANN_CONSTANT*temperature*(number_of_links as f64)*parameters.abs_tol);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&force, &temperature);
            let gibbs_free_energy_per_link_0 = model.gibbs_free_energy_per_link(&0.0, &temperature);
            let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&force, &temperature);
            let residual_abs = &gibbs_free_energy_per_link - &gibbs_free_energy_per_link_0 - &relative_gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&gibbs_free_energy_per_link;
            assert!(residual_abs.abs() <= BOLTZMANN_CONSTANT*temperature*parameters.abs_tol);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature);
            let nondimensional_gibbs_free_energy_0 = model.nondimensional_gibbs_free_energy(&0.0, &temperature);
            let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
            let residual_abs = &nondimensional_gibbs_free_energy - &nondimensional_gibbs_free_energy_0 - &nondimensional_relative_gibbs_free_energy;
            let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy;
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature);
            let nondimensional_gibbs_free_energy_per_link_0 = model.nondimensional_gibbs_free_energy_per_link(&0.0, &temperature);
            let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
            let residual_abs = &nondimensional_gibbs_free_energy_per_link - &nondimensional_gibbs_free_energy_per_link_0 - &nondimensional_relative_gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link;
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
    fn end_to_end_length()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length_0 = model.end_to_end_length(&0.0, &temperature);
            assert!(end_to_end_length_0.abs() <= 3.0*BOLTZMANN_CONSTANT*temperature*(number_of_links as f64)*0.0);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length_per_link_0 = model.end_to_end_length_per_link(&0.0, &temperature);
            assert!(end_to_end_length_per_link_0.abs() <= 3.0*BOLTZMANN_CONSTANT*temperature*0.0);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_0 = model.nondimensional_end_to_end_length(&0.0);
            assert!(nondimensional_end_to_end_length_0.abs() <= 3.0*(number_of_links as f64)*0.0);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link_0 = model.nondimensional_end_to_end_length_per_link(&0.0);
            assert!(nondimensional_end_to_end_length_per_link_0.abs() <= 3.0*0.0);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_gibbs_free_energy_0 = model.relative_gibbs_free_energy(&0.0, &temperature);
            assert!(relative_gibbs_free_energy_0.abs() <= 0.0);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_gibbs_free_energy_per_link_0 = model.relative_gibbs_free_energy_per_link(&0.0, &temperature);
            assert!(relative_gibbs_free_energy_per_link_0.abs() <= 0.0);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_relative_gibbs_free_energy_0 = model.nondimensional_relative_gibbs_free_energy(&0.0);
            assert!(nondimensional_relative_gibbs_free_energy_0.abs() <= 0.0);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_relative_gibbs_free_energy_per_link_0 = model.nondimensional_relative_gibbs_free_energy_per_link(&0.0);
            assert!(nondimensional_relative_gibbs_free_energy_per_link_0.abs() <= 0.0);
        }
    }
}
mod connection
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + 0.5*parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.end_to_end_length(&force, &temperature);
            let h = parameters.rel_tol*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length_from_derivative = -(model.relative_gibbs_free_energy(&(force + 0.5*h), &temperature) - model.relative_gibbs_free_energy(&(force - 0.5*h), &temperature))/h;
            let residual_abs = &end_to_end_length - &end_to_end_length_from_derivative;
            let residual_rel = &residual_abs/&end_to_end_length;
            assert!(residual_rel.abs() <= h);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + 0.5*parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length_per_link = model.end_to_end_length_per_link(&force, &temperature);
            let h = parameters.rel_tol*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length_per_link_from_derivative = -(model.relative_gibbs_free_energy_per_link(&(force + 0.5*h), &temperature) - model.relative_gibbs_free_energy_per_link(&(force - 0.5*h), &temperature))/h;
            let residual_abs = &end_to_end_length_per_link - &end_to_end_length_per_link_from_derivative;
            let residual_rel = &residual_abs/&end_to_end_length_per_link;
            assert!(residual_rel.abs() <= h);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + 0.5*parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_force);
            let h = parameters.rel_tol;
            let nondimensional_end_to_end_length_from_derivative = -(model.nondimensional_relative_gibbs_free_energy(&(nondimensional_force + 0.5*h)) - model.nondimensional_relative_gibbs_free_energy(&(nondimensional_force - 0.5*h)))/h;
            let residual_abs = &nondimensional_end_to_end_length - &nondimensional_end_to_end_length_from_derivative;
            let residual_rel = &residual_abs/&nondimensional_end_to_end_length;
            assert!(residual_rel.abs() <= h);
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
            let model = Ideal::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + 0.5*parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&nondimensional_force);
            let h = parameters.rel_tol;
            let nondimensional_end_to_end_length_per_link_from_derivative = -(model.nondimensional_relative_gibbs_free_energy_per_link(&(nondimensional_force + 0.5*h)) - model.nondimensional_relative_gibbs_free_energy_per_link(&(nondimensional_force - 0.5*h)))/h;
            let residual_abs = &nondimensional_end_to_end_length_per_link - &nondimensional_end_to_end_length_per_link_from_derivative;
            let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link;
            assert!(residual_rel.abs() <= h);
        }
    }
}
