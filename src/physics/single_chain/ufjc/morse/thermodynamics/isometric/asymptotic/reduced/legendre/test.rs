#![cfg(test)]
use super::*;
use crate::physics::single_chain::test::Parameters;
use crate::physics::single_chain::ufjc::morse::thermodynamics::isotensional::asymptotic::reduced::nondimensional_end_to_end_length_per_link as isotensional_nondimensional_end_to_end_length_per_link;
mod base
{
    use super::*;
    use rand::Rng;
    #[test]
    fn init()
    {
        let parameters = Parameters::default();
        let _ = MORSEFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference, parameters.link_energy_reference);
    }
    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            assert_eq!(number_of_links, MORSEFJC::init(number_of_links, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference, parameters.link_energy_reference).number_of_links);
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
            assert_eq!(link_length, MORSEFJC::init(parameters.number_of_links_minimum, link_length, parameters.hinge_mass_reference, parameters.link_stiffness_reference, parameters.link_energy_reference).link_length);
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
            assert_eq!(hinge_mass, MORSEFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, hinge_mass, parameters.link_stiffness_reference, parameters.link_energy_reference).hinge_mass);
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
            assert_eq!(link_stiffness, MORSEFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, link_stiffness, parameters.link_energy_reference).link_stiffness);
        }
    }
    #[test]
    fn link_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(link_energy, MORSEFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference, link_energy).link_energy);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            assert_eq!(number_of_links, model.number_of_links);
            assert_eq!(link_length, model.link_length);
            assert_eq!(hinge_mass, model.hinge_mass);
            assert_eq!(link_stiffness, model.link_stiffness);
            assert_eq!(link_energy, model.link_energy);
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
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let nondimensional_force = model.nondimensional_force(&nondimensional_end_to_end_length_per_link, &temperature);
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let force = model.force(&end_to_end_length, &temperature);
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
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let helmholtz_free_energy = model.helmholtz_free_energy(&end_to_end_length, &temperature);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&end_to_end_length, &temperature);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let helmholtz_free_energy = model.helmholtz_free_energy(&end_to_end_length, &temperature);
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&end_to_end_length, &temperature);
            let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let helmholtz_free_energy = model.helmholtz_free_energy(&end_to_end_length, &temperature);
            let helmholtz_free_energy_0 = model.helmholtz_free_energy(&(ZERO*(number_of_links as f64)*link_length), &temperature);
            let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&end_to_end_length, &temperature);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
            let helmholtz_free_energy_per_link_0 = model.helmholtz_free_energy_per_link(&(ZERO*(number_of_links as f64)*link_length), &temperature);
            let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let nondimensional_helmholtz_free_energy_0 = model.nondimensional_helmholtz_free_energy(&ZERO, &temperature);
            let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
            let nondimensional_helmholtz_free_energy_per_link_0 = model.nondimensional_helmholtz_free_energy_per_link(&ZERO, &temperature);
            let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
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
    use crate::physics::single_chain::ZERO;
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
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force_0 = model.force(&(ZERO*(number_of_links as f64)*link_length), &temperature);
            assert!(force_0.abs() <= 3.1*BOLTZMANN_CONSTANT*temperature/link_length*ZERO);
        }
    }
    #[test]
    fn nondimensional_force()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_0 = model.nondimensional_force(&ZERO, &temperature);
            assert!(nondimensional_force_0.abs() <= 3.1*ZERO);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_helmholtz_free_energy_0 = model.relative_helmholtz_free_energy(&(ZERO*(number_of_links as f64)*link_length), &temperature);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_helmholtz_free_energy_per_link_0 = model.relative_helmholtz_free_energy_per_link(&(ZERO*(number_of_links as f64)*link_length), &temperature);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
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
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_helmholtz_free_energy_per_link_0 = model.nondimensional_relative_helmholtz_free_energy_per_link(&ZERO, &temperature);
            assert!(nondimensional_relative_helmholtz_free_energy_per_link_0.abs() <= ZERO);
        }
    }
}
mod connection
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
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let force = model.force(&end_to_end_length, &temperature);
            let h = parameters.rel_tol*(number_of_links as f64)*link_length;
            let force_from_derivative = (model.relative_helmholtz_free_energy(&(end_to_end_length + 0.5*h), &temperature) - model.relative_helmholtz_free_energy(&(end_to_end_length - 0.5*h), &temperature))/h;
            let residual_abs = &force - &force_from_derivative;
            let residual_rel = &residual_abs/&force;
            assert!(residual_rel.abs() <= h);
        }
    }
    #[test]
    fn nondimensional_force()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let link_energy = parameters.link_energy_reference + parameters.link_energy_scale*(0.5 - rng.gen::<f64>());
            let model = MORSEFJC::init(number_of_links, link_length, hinge_mass, link_stiffness, link_energy);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force_max = (link_stiffness*link_energy/8.0).sqrt()/BOLTZMANN_CONSTANT/temperature*link_length;
            let nondimensional_end_to_end_length_per_link_max = isotensional_nondimensional_end_to_end_length_per_link( &(link_stiffness*link_length.powi(2)/BOLTZMANN_CONSTANT/temperature), &(link_energy/BOLTZMANN_CONSTANT/temperature), &(0.999*nondimensional_force_max));
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length_per_link_max*rng.gen::<f64>();
            let nondimensional_force = model.nondimensional_force(&nondimensional_end_to_end_length_per_link, &temperature);
            let h = parameters.rel_tol;
            let nondimensional_force_from_derivative = (model.nondimensional_relative_helmholtz_free_energy_per_link(&(nondimensional_end_to_end_length_per_link + 0.5*h), &temperature) - model.nondimensional_relative_helmholtz_free_energy_per_link(&(nondimensional_end_to_end_length_per_link - 0.5*h), &temperature))/h;
            let residual_abs = &nondimensional_force - &nondimensional_force_from_derivative;
            let residual_rel = &residual_abs/&nondimensional_force;
            assert!(residual_rel.abs() <= h);
        }
    }
}
