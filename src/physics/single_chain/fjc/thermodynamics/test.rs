#![cfg(test)]
use super::*;
use crate::physics::BOLTZMANN_CONSTANT;
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
mod thermodynamic_limit
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
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let force = model.isometric.force(&end_to_end_length, &temperature);
            let end_to_end_length_out = model.isotensional.end_to_end_length(&force, &temperature);
            let residual_abs = &end_to_end_length - &end_to_end_length_out;
            let residual_rel = &residual_abs/&end_to_end_length;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn end_to_end_length_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let force = model.isometric.force(&end_to_end_length, &temperature);
            let end_to_end_length_per_link = nondimensional_end_to_end_length_per_link*link_length;
            let end_to_end_length_per_link_out = model.isotensional.end_to_end_length_per_link(&force, &temperature);
            let residual_abs = &end_to_end_length_per_link - &end_to_end_length_per_link_out;
            let residual_rel = &residual_abs/&end_to_end_length_per_link;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn nondimensional_end_to_end_length()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force = model.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link);
            let nondimensional_end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64);
            let nondimensional_end_to_end_length_out = model.isotensional.nondimensional_end_to_end_length(&nondimensional_force);
            let residual_abs = &nondimensional_end_to_end_length - &nondimensional_end_to_end_length_out;
            let residual_rel = &residual_abs/&nondimensional_end_to_end_length;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn nondimensional_end_to_end_length_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force = model.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link);
            let nondimensional_end_to_end_length_per_link_out = model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
            let residual_abs = &nondimensional_end_to_end_length_per_link - &nondimensional_end_to_end_length_per_link_out;
            let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn force()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.end_to_end_length(&force, &temperature);
            let force_out = model.isometric.force(&end_to_end_length, &temperature);
            let residual_abs = &force - &force_out;
            let residual_rel = &residual_abs/&force;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn nondimensional_force()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length_per_link = model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
            let nondimensional_force_out = model.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link);
            let residual_abs = &nondimensional_force - &nondimensional_force_out;
            let residual_rel = &residual_abs/&nondimensional_force;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let force = model.isometric.force(&end_to_end_length, &temperature);
            let helmholtz_free_energy = model.isometric.helmholtz_free_energy(&end_to_end_length, &temperature);
            let helmholtz_free_energy_out = model.isotensional.gibbs_free_energy(&force, &temperature) + force*end_to_end_length;
            let residual_abs = &helmholtz_free_energy - &helmholtz_free_energy_out;
            let residual_rel = &residual_abs/&helmholtz_free_energy;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn helmholtz_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let end_to_end_length_per_link = nondimensional_end_to_end_length_per_link*link_length;
            let force = model.isometric.force(&end_to_end_length, &temperature);
            let helmholtz_free_energy_per_link = model.isometric.helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
            let helmholtz_free_energy_per_link_out = model.isotensional.gibbs_free_energy_per_link(&force, &temperature) + force*end_to_end_length_per_link;
            let residual_abs = &helmholtz_free_energy_per_link - &helmholtz_free_energy_per_link_out;
            let residual_rel = &residual_abs/&helmholtz_free_energy_per_link;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn relative_helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let force = model.isometric.force(&end_to_end_length, &temperature);
            let relative_helmholtz_free_energy = model.isometric.relative_helmholtz_free_energy(&end_to_end_length, &temperature);
            let relative_helmholtz_free_energy_out = model.isotensional.relative_gibbs_free_energy(&force, &temperature) + force*end_to_end_length;
            let residual_abs = &relative_helmholtz_free_energy - &relative_helmholtz_free_energy_out;
            let residual_rel = &residual_abs/&relative_helmholtz_free_energy;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn relative_helmholtz_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let end_to_end_length_per_link = nondimensional_end_to_end_length_per_link*link_length;
            let force = model.isometric.force(&end_to_end_length, &temperature);
            let relative_helmholtz_free_energy_per_link = model.isometric.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
            let relative_helmholtz_free_energy_per_link_out = model.isotensional.relative_gibbs_free_energy_per_link(&force, &temperature) + force*end_to_end_length_per_link;
            let residual_abs = &relative_helmholtz_free_energy_per_link - &relative_helmholtz_free_energy_per_link_out;
            let residual_rel = &residual_abs/&relative_helmholtz_free_energy_per_link;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn nondimensional_helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64);
            let nondimensional_force = model.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link);
            let nondimensional_helmholtz_free_energy = model.isometric.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let nondimensional_helmholtz_free_energy_out = model.isotensional.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature) + nondimensional_force*nondimensional_end_to_end_length;
            let residual_abs = &nondimensional_helmholtz_free_energy - &nondimensional_helmholtz_free_energy_out;
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn nondimensional_helmholtz_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force = model.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link);
            let nondimensional_helmholtz_free_energy_per_link = model.isometric.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
            let nondimensional_helmholtz_free_energy_per_link_out = model.isotensional.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature) + nondimensional_force*nondimensional_end_to_end_length_per_link;
            let residual_abs = &nondimensional_helmholtz_free_energy_per_link - &nondimensional_helmholtz_free_energy_per_link_out;
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn nondimensional_relative_helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64);
            let nondimensional_force = model.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link);
            let nondimensional_relative_helmholtz_free_energy = model.isometric.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link);
            let nondimensional_relative_helmholtz_free_energy_out = model.isotensional.nondimensional_relative_gibbs_free_energy(&nondimensional_force) + nondimensional_force*nondimensional_end_to_end_length;
            let residual_abs = &nondimensional_relative_helmholtz_free_energy - &nondimensional_relative_helmholtz_free_energy_out;
            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn nondimensional_relative_helmholtz_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force = model.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link);
            let nondimensional_relative_helmholtz_free_energy_per_link = model.isometric.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
            let nondimensional_relative_helmholtz_free_energy_per_link_out = model.isotensional.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force) + nondimensional_force*nondimensional_end_to_end_length_per_link;
            let residual_abs = &nondimensional_relative_helmholtz_free_energy_per_link - &nondimensional_relative_helmholtz_free_energy_per_link_out;
            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.end_to_end_length(&force, &temperature);
            let gibbs_free_energy = model.isotensional.gibbs_free_energy(&force, &temperature);
            let gibbs_free_energy_out = model.isometric.helmholtz_free_energy(&end_to_end_length, &temperature) - force*end_to_end_length;
            let residual_abs = &gibbs_free_energy - &gibbs_free_energy_out;
            let residual_rel = &residual_abs/&gibbs_free_energy;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.end_to_end_length(&force, &temperature);
            let end_to_end_length_per_link = end_to_end_length/(number_of_links as f64);
            let gibbs_free_energy_per_link = model.isotensional.gibbs_free_energy_per_link(&force, &temperature);
            let gibbs_free_energy_per_link_out = model.isometric.helmholtz_free_energy_per_link(&end_to_end_length, &temperature) - force*end_to_end_length_per_link;
            let residual_abs = &gibbs_free_energy_per_link - &gibbs_free_energy_per_link_out;
            let residual_rel = &residual_abs/&gibbs_free_energy_per_link;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn relative_gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.end_to_end_length(&force, &temperature);
            let relative_gibbs_free_energy = model.isotensional.relative_gibbs_free_energy(&force, &temperature);
            let relative_gibbs_free_energy_out = model.isometric.relative_helmholtz_free_energy(&end_to_end_length, &temperature) - force*end_to_end_length;
            let residual_abs = &relative_gibbs_free_energy - &relative_gibbs_free_energy_out;
            let residual_rel = &residual_abs/&relative_gibbs_free_energy;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn relative_gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.end_to_end_length(&force, &temperature);
            let end_to_end_length_per_link = end_to_end_length/(number_of_links as f64);
            let relative_gibbs_free_energy_per_link = model.isotensional.relative_gibbs_free_energy_per_link(&force, &temperature);
            let relative_gibbs_free_energy_per_link_out = model.isometric.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature) - force*end_to_end_length_per_link;
            let residual_abs = &relative_gibbs_free_energy_per_link - &relative_gibbs_free_energy_per_link_out;
            let residual_rel = &residual_abs/&relative_gibbs_free_energy_per_link;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn nondimensional_gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = model.isotensional.nondimensional_end_to_end_length(&nondimensional_force);
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length/(number_of_links as f64);
            let nondimensional_gibbs_free_energy = model.isotensional.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature);
            let nondimensional_gibbs_free_energy_out = model.isometric.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature) - nondimensional_force*nondimensional_end_to_end_length;
            let residual_abs = &nondimensional_gibbs_free_energy - &nondimensional_gibbs_free_energy_out;
            let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn nondimensional_gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length_per_link = model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
            let nondimensional_gibbs_free_energy_per_link = model.isotensional.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature);
            let nondimensional_gibbs_free_energy_per_link_out = model.isometric.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature) - nondimensional_force*nondimensional_end_to_end_length_per_link;
            let residual_abs = &nondimensional_gibbs_free_energy_per_link - &nondimensional_gibbs_free_energy_per_link_out;
            let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn nondimensional_relative_gibbs_free_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = model.isotensional.nondimensional_end_to_end_length(&nondimensional_force);
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length/(number_of_links as f64);
            let nondimensional_relative_gibbs_free_energy = model.isotensional.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
            let nondimensional_relative_gibbs_free_energy_out = model.isometric.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link) - nondimensional_force*nondimensional_end_to_end_length;
            let residual_abs = &nondimensional_relative_gibbs_free_energy - &nondimensional_relative_gibbs_free_energy_out;
            let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
    #[test]
    fn nondimensional_relative_gibbs_free_energy_per_link()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = parameters.number_of_links_maximum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length_per_link = model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
            let nondimensional_relative_gibbs_free_energy_per_link = model.isotensional.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
            let nondimensional_relative_gibbs_free_energy_per_link_out = model.isometric.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link) - nondimensional_force*nondimensional_end_to_end_length_per_link;
            let residual_abs = &nondimensional_relative_gibbs_free_energy_per_link - &nondimensional_relative_gibbs_free_energy_per_link_out;
            let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
            assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
        }
    }
}
mod modified_canonical
{
    use super::*;
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
                    let integrand_numerator = |end_to_end_length: f64|
                    {
                        (model.modified_canonical.force(&end_to_end_length, &potential_stiffness, &temperature) - model.isometric.force(&end_to_end_length, &temperature)).powi(2)
                    };
                    let integrand_denominator = |end_to_end_length: f64|
                    {
                        model.modified_canonical.force(&end_to_end_length, &potential_stiffness, &temperature).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                    let denominator = integrate(integrand_denominator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                    let integrand_numerator = |nondimensional_end_to_end_length_per_link: f64|
                    {
                        (model.modified_canonical.nondimensional_force(&nondimensional_end_to_end_length_per_link, &nondimensional_potential_stiffness) - model.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_end_to_end_length_per_link: f64|
                    {
                        model.modified_canonical.nondimensional_force(&nondimensional_end_to_end_length_per_link, &nondimensional_potential_stiffness).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                    let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                    let integrand_numerator = |end_to_end_length: f64|
                    {
                        (model.modified_canonical.relative_helmholtz_free_energy(&end_to_end_length, &potential_stiffness, &temperature) - model.isometric.relative_helmholtz_free_energy(&end_to_end_length, &temperature)).powi(2)
                    };
                    let integrand_denominator = |end_to_end_length: f64|
                    {
                        model.modified_canonical.relative_helmholtz_free_energy(&end_to_end_length, &potential_stiffness, &temperature).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                    let denominator = integrate(integrand_denominator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                    let integrand_numerator = |end_to_end_length: f64|
                    {
                        (model.modified_canonical.relative_helmholtz_free_energy_per_link(&end_to_end_length, &potential_stiffness, &temperature) - model.isometric.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature)).powi(2)
                    };
                    let integrand_denominator = |end_to_end_length: f64|
                    {
                        model.modified_canonical.relative_helmholtz_free_energy_per_link(&end_to_end_length, &potential_stiffness, &temperature).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                    let denominator = integrate(integrand_denominator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                        (model.modified_canonical.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &nondimensional_potential_stiffness) - model.isometric.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_end_to_end_length_per_link: f64|
                    {
                        model.modified_canonical.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &nondimensional_potential_stiffness).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                    let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                        (model.modified_canonical.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &nondimensional_potential_stiffness) - model.isometric.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_end_to_end_length_per_link: f64|
                    {
                        model.modified_canonical.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &nondimensional_potential_stiffness).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                    let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
            }
        }
    }
    mod weak_potential
    {
        use super::*;
        use rand::Rng;
        use crate::physics::single_chain::test::integrate;
        use crate::physics::single_chain::POINTS;
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
                        let force = model.modified_canonical.force(&potential_distance, &potential_stiffness, &temperature);
                        (model.modified_canonical.end_to_end_length(&potential_distance, &potential_stiffness, &temperature) - model.isotensional.end_to_end_length(&force, &temperature)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_potential_distance: f64|
                    {
                        let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                        model.modified_canonical.end_to_end_length(&potential_distance, &potential_stiffness, &temperature).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                        let force = model.modified_canonical.force(&potential_distance, &potential_stiffness, &temperature);
                        (model.modified_canonical.end_to_end_length_per_link(&potential_distance, &potential_stiffness, &temperature) - model.isotensional.end_to_end_length_per_link(&force, &temperature)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_potential_distance: f64|
                    {
                        let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                        model.modified_canonical.end_to_end_length_per_link(&potential_distance, &potential_stiffness, &temperature).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                        let nondimensional_force = model.modified_canonical.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                        (model.modified_canonical.nondimensional_end_to_end_length(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.isotensional.nondimensional_end_to_end_length(&nondimensional_force)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_potential_distance: f64|
                    {
                        model.modified_canonical.nondimensional_end_to_end_length(&nondimensional_potential_distance, &nondimensional_potential_stiffness).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                        let nondimensional_force = model.modified_canonical.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                        (model.modified_canonical.nondimensional_end_to_end_length_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_potential_distance: f64|
                    {
                        model.modified_canonical.nondimensional_end_to_end_length_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                let potential_distance_ref = parameters.nondimensional_potential_distance_large_1*(number_of_links as f64)*link_length;
                let residual_rel = |nondimensional_potential_stiffness|
                {
                    let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                    let force_ref = model.modified_canonical.force(&potential_distance_ref, &potential_stiffness, &temperature);
                    let integrand_numerator = |nondimensional_potential_distance: f64|
                    {
                        let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                        let force = model.modified_canonical.force(&potential_distance, &potential_stiffness, &temperature);
                        (model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature) - model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy(&potential_distance_ref, &potential_stiffness, &temperature) - model.isotensional.gibbs_free_energy(&force, &temperature) + model.isotensional.gibbs_free_energy(&force_ref, &temperature)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_potential_distance: f64|
                    {
                        let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                        (model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature) - model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy(&potential_distance_ref, &potential_stiffness, &temperature)).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                let potential_distance_ref = parameters.nondimensional_potential_distance_large_1*(number_of_links as f64)*link_length;
                let residual_rel = |nondimensional_potential_stiffness|
                {
                    let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                    let force_ref = model.modified_canonical.force(&potential_distance_ref, &potential_stiffness, &temperature);
                    let integrand_numerator = |nondimensional_potential_distance: f64|
                    {
                        let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                        let force = model.modified_canonical.force(&potential_distance, &potential_stiffness, &temperature);
                        (model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) - model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy_per_link(&potential_distance_ref, &potential_stiffness, &temperature) - model.isotensional.gibbs_free_energy_per_link(&force, &temperature) + model.isotensional.gibbs_free_energy_per_link(&force_ref, &temperature)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_potential_distance: f64|
                    {
                        let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                        (model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) - model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy_per_link(&potential_distance_ref, &potential_stiffness, &temperature)).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                let potential_distance_ref = parameters.nondimensional_potential_distance_large_1*(number_of_links as f64)*link_length;
                let residual_rel = |nondimensional_potential_stiffness|
                {
                    let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                    let force_ref = model.modified_canonical.force(&potential_distance_ref, &potential_stiffness, &temperature);
                    let integrand_numerator = |nondimensional_potential_distance: f64|
                    {
                        let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                        let force = model.modified_canonical.force(&potential_distance, &potential_stiffness, &temperature);
                        (model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature) - model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy(&potential_distance_ref, &potential_stiffness, &temperature) - model.isotensional.relative_gibbs_free_energy(&force, &temperature) + model.isotensional.relative_gibbs_free_energy(&force_ref, &temperature)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_potential_distance: f64|
                    {
                        let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                        (model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature) - model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy(&potential_distance_ref, &potential_stiffness, &temperature)).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                let potential_distance_ref = parameters.nondimensional_potential_distance_large_1*(number_of_links as f64)*link_length;
                let residual_rel = |nondimensional_potential_stiffness|
                {
                    let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                    let force_ref = model.modified_canonical.force(&potential_distance_ref, &potential_stiffness, &temperature);
                    let integrand_numerator = |nondimensional_potential_distance: f64|
                    {
                        let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                        let force = model.modified_canonical.force(&potential_distance, &potential_stiffness, &temperature);
                        (model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) - model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy_per_link(&potential_distance_ref, &potential_stiffness, &temperature) - model.isotensional.relative_gibbs_free_energy_per_link(&force, &temperature) + model.isotensional.relative_gibbs_free_energy_per_link(&force_ref, &temperature)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_potential_distance: f64|
                    {
                        let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                        (model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) - model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy_per_link(&potential_distance_ref, &potential_stiffness, &temperature)).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                let nondimensional_potential_distance_ref = parameters.nondimensional_potential_distance_large_1;
                let residual_rel = |nondimensional_potential_stiffness|
                {
                    let nondimensional_force_ref = model.modified_canonical.nondimensional_force(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness);
                    let integrand_numerator = |nondimensional_potential_distance: f64|
                    {
                        let nondimensional_force = model.modified_canonical.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                        (model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) - model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness, &temperature) - model.isotensional.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature) + model.isotensional.nondimensional_gibbs_free_energy(&nondimensional_force_ref, &temperature)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_potential_distance: f64|
                    {
                        (model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) - model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness, &temperature)).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                let nondimensional_potential_distance_ref = parameters.nondimensional_potential_distance_large_1;
                let residual_rel = |nondimensional_potential_stiffness|
                {
                    let nondimensional_force_ref = model.modified_canonical.nondimensional_force(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness);
                    let integrand_numerator = |nondimensional_potential_distance: f64|
                    {
                        let nondimensional_force = model.modified_canonical.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                        (model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) - model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness, &temperature) - model.isotensional.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature) + model.isotensional.nondimensional_gibbs_free_energy_per_link(&nondimensional_force_ref, &temperature)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_potential_distance: f64|
                    {
                        (model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) - model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness, &temperature)).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                let nondimensional_potential_distance_ref = parameters.nondimensional_potential_distance_large_1;
                let residual_rel = |nondimensional_potential_stiffness|
                {
                    let nondimensional_force_ref = model.modified_canonical.nondimensional_force(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness);
                    let integrand_numerator = |nondimensional_potential_distance: f64|
                    {
                        let nondimensional_force = model.modified_canonical.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                        (model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness) - model.isotensional.nondimensional_relative_gibbs_free_energy(&nondimensional_force) + model.isotensional.nondimensional_relative_gibbs_free_energy(&nondimensional_force_ref)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_potential_distance: f64|
                    {
                        (model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness)).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                let nondimensional_potential_distance_ref = parameters.nondimensional_potential_distance_large_1;
                let residual_rel = |nondimensional_potential_stiffness|
                {
                    let nondimensional_force_ref = model.modified_canonical.nondimensional_force(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness);
                    let integrand_numerator = |nondimensional_potential_distance: f64|
                    {
                        let nondimensional_force = model.modified_canonical.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                        (model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness) - model.isotensional.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force) + model.isotensional.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force_ref)).powi(2)
                    };
                    let integrand_denominator = |nondimensional_potential_distance: f64|
                    {
                        (model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness)).powi(2)
                    };
                    let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                    (numerator/denominator).sqrt()
                };
                let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
            }
        }
    }
    mod asymptotic
    {
        use super::*;
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
                        let integrand_numerator = |end_to_end_length: f64|
                        {
                            (model.isometric.force(&end_to_end_length, &temperature) -model.modified_canonical.asymptotic.strong_potential.force(&end_to_end_length, &potential_stiffness, &temperature)).powi(2)
                        };
                        let integrand_denominator = |end_to_end_length: f64|
                        {
                            model.isometric.force(&end_to_end_length, &temperature).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                        let denominator = integrate(integrand_denominator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                    let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                        let integrand_numerator = |nondimensional_end_to_end_length_per_link: f64|
                        {
                            (model.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link) - model.modified_canonical.asymptotic.strong_potential.nondimensional_force(&nondimensional_end_to_end_length_per_link, &nondimensional_potential_stiffness)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_end_to_end_length_per_link: f64|
                        {
                            model.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                        let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                    let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                        let integrand_numerator = |end_to_end_length: f64|
                        {
                            (model.isometric.helmholtz_free_energy(&end_to_end_length, &temperature) - model.isometric.helmholtz_free_energy(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &temperature) - model.modified_canonical.asymptotic.strong_potential.helmholtz_free_energy(&end_to_end_length, &potential_stiffness, &temperature) + model.modified_canonical.asymptotic.strong_potential.helmholtz_free_energy(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature)).powi(2)
                        };
                        let integrand_denominator = |end_to_end_length: f64|
                        {
                            (model.isometric.helmholtz_free_energy(&end_to_end_length, &temperature) - model.isometric.helmholtz_free_energy(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &temperature)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                        let denominator = integrate(integrand_denominator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                    let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                        let integrand_numerator = |end_to_end_length: f64|
                        {
                            (model.isometric.helmholtz_free_energy_per_link(&end_to_end_length, &temperature) - model.isometric.helmholtz_free_energy_per_link(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &temperature) - model.modified_canonical.asymptotic.strong_potential.helmholtz_free_energy_per_link(&end_to_end_length, &potential_stiffness, &temperature) + model.modified_canonical.asymptotic.strong_potential.helmholtz_free_energy_per_link(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature)).powi(2)
                        };
                        let integrand_denominator = |end_to_end_length: f64|
                        {
                            (model.isometric.helmholtz_free_energy_per_link(&end_to_end_length, &temperature) - model.isometric.helmholtz_free_energy_per_link(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &temperature)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                        let denominator = integrate(integrand_denominator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                    let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                        let integrand_numerator = |end_to_end_length: f64|
                        {
                            (model.isometric.relative_helmholtz_free_energy(&end_to_end_length, &temperature) - model.isometric.relative_helmholtz_free_energy(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &temperature) - model.modified_canonical.asymptotic.strong_potential.relative_helmholtz_free_energy(&end_to_end_length, &potential_stiffness, &temperature) + model.modified_canonical.asymptotic.strong_potential.relative_helmholtz_free_energy(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature)).powi(2)
                        };
                        let integrand_denominator = |end_to_end_length: f64|
                        {
                            (model.isometric.relative_helmholtz_free_energy(&end_to_end_length, &temperature) - model.isometric.relative_helmholtz_free_energy(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &temperature)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                        let denominator = integrate(integrand_denominator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                    let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                        let integrand_numerator = |end_to_end_length: f64|
                        {
                            (model.isometric.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature) - model.isometric.relative_helmholtz_free_energy_per_link(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &temperature) - model.modified_canonical.asymptotic.strong_potential.relative_helmholtz_free_energy_per_link(&end_to_end_length, &potential_stiffness, &temperature) + model.modified_canonical.asymptotic.strong_potential.relative_helmholtz_free_energy_per_link(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &potential_stiffness, &temperature)).powi(2)
                        };
                        let integrand_denominator = |end_to_end_length: f64|
                        {
                            (model.isometric.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature) - model.isometric.relative_helmholtz_free_energy_per_link(&(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &temperature)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                        let denominator = integrate(integrand_denominator, &(ZERO*(number_of_links as f64)*link_length), &(parameters.nondimensional_potential_distance_small*(number_of_links as f64)*link_length), &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                    let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                        let integrand_numerator = |nondimensional_end_to_end_length_per_link: f64|
                        {
                            (model.isometric.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature) - model.isometric.nondimensional_helmholtz_free_energy(&parameters.nondimensional_potential_distance_small, &temperature) - model.modified_canonical.asymptotic.strong_potential.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &nondimensional_potential_stiffness, &temperature) + model.modified_canonical.asymptotic.strong_potential.nondimensional_helmholtz_free_energy(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness, &temperature)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_end_to_end_length_per_link: f64|
                        {
                            (model.isometric.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature) - model.isometric.nondimensional_helmholtz_free_energy(&parameters.nondimensional_potential_distance_small, &temperature)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                        let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                    let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                        let integrand_numerator = |nondimensional_end_to_end_length_per_link: f64|
                        {
                            (model.isometric.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature) - model.isometric.nondimensional_helmholtz_free_energy_per_link(&parameters.nondimensional_potential_distance_small, &temperature) - model.modified_canonical.asymptotic.strong_potential.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &nondimensional_potential_stiffness, &temperature) + model.modified_canonical.asymptotic.strong_potential.nondimensional_helmholtz_free_energy_per_link(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness, &temperature)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_end_to_end_length_per_link: f64|
                        {
                            (model.isometric.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature) - model.isometric.nondimensional_helmholtz_free_energy_per_link(&parameters.nondimensional_potential_distance_small, &temperature)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                        let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                    let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                            (model.isometric.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link) - model.isometric.nondimensional_relative_helmholtz_free_energy(&parameters.nondimensional_potential_distance_small) - model.modified_canonical.asymptotic.strong_potential.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &nondimensional_potential_stiffness) + model.modified_canonical.asymptotic.strong_potential.nondimensional_relative_helmholtz_free_energy(&parameters.nondimensional_potential_distance_small, &nondimensional_potential_stiffness)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_end_to_end_length_per_link: f64|
                        {
                            (model.isometric.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link) - model.isometric.nondimensional_relative_helmholtz_free_energy(&parameters.nondimensional_potential_distance_small)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                        let denominator = integrate(integrand_denominator, &ZERO, &parameters.nondimensional_potential_distance_small, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_large);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_large*parameters.log_log_scale);
                    let log_log_slope = (residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
                }
            }
        }
        mod weak_potential
        {
            use super::*;
            use rand::Rng;
            use crate::physics::single_chain::test::integrate;
            use crate::physics::single_chain::POINTS;
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
                            let force = model.modified_canonical.asymptotic.weak_potential.force(&potential_distance, &potential_stiffness);
                            (model.modified_canonical.asymptotic.weak_potential.end_to_end_length(&potential_distance, &potential_stiffness, &temperature) - model.isotensional.end_to_end_length(&force, &temperature)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_potential_distance: f64|
                        {
                            let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                            model.modified_canonical.asymptotic.weak_potential.end_to_end_length(&potential_distance, &potential_stiffness, &temperature).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                    let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                    assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                            let force = model.modified_canonical.asymptotic.weak_potential.force(&potential_distance, &potential_stiffness);
                            (model.modified_canonical.asymptotic.weak_potential.end_to_end_length_per_link(&potential_distance, &potential_stiffness, &temperature) - model.isotensional.end_to_end_length_per_link(&force, &temperature)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_potential_distance: f64|
                        {
                            let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                            model.modified_canonical.asymptotic.weak_potential.end_to_end_length_per_link(&potential_distance, &potential_stiffness, &temperature).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                    let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                    assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                            let nondimensional_force = model.modified_canonical.asymptotic.weak_potential.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                            (model.modified_canonical.asymptotic.weak_potential.nondimensional_end_to_end_length(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.isotensional.nondimensional_end_to_end_length(&nondimensional_force)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_potential_distance: f64|
                        {
                            model.modified_canonical.asymptotic.weak_potential.nondimensional_end_to_end_length(&nondimensional_potential_distance, &nondimensional_potential_stiffness).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                    let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                    assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                            let nondimensional_force = model.modified_canonical.asymptotic.weak_potential.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                            (model.modified_canonical.asymptotic.weak_potential.nondimensional_end_to_end_length_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness) - model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_potential_distance: f64|
                        {
                            model.modified_canonical.asymptotic.weak_potential.nondimensional_end_to_end_length_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                    let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                    assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                    let potential_distance_ref = parameters.nondimensional_potential_distance_large_1*(number_of_links as f64)*link_length;
                    let residual_rel = |nondimensional_potential_stiffness|
                    {
                        let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                        let force_ref = model.modified_canonical.asymptotic.weak_potential.force(&potential_distance_ref, &potential_stiffness);
                        let integrand_numerator = |nondimensional_potential_distance: f64|
                        {
                            let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                            let force = model.modified_canonical.asymptotic.weak_potential.force(&potential_distance, &potential_stiffness);
                            (model.isotensional.gibbs_free_energy(&force, &temperature) - model.isotensional.gibbs_free_energy(&force_ref, &temperature) - model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature) + model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy(&potential_distance_ref, &potential_stiffness, &temperature)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_potential_distance: f64|
                        {
                            let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                            let force = model.modified_canonical.asymptotic.weak_potential.force(&potential_distance, &potential_stiffness);
                            (model.isotensional.gibbs_free_energy(&force, &temperature) - model.isotensional.gibbs_free_energy(&force_ref, &temperature)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                    let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                    assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                    let potential_distance_ref = parameters.nondimensional_potential_distance_large_1*(number_of_links as f64)*link_length;
                    let residual_rel = |nondimensional_potential_stiffness|
                    {
                        let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                        let force_ref = model.modified_canonical.asymptotic.weak_potential.force(&potential_distance_ref, &potential_stiffness);
                        let integrand_numerator = |nondimensional_potential_distance: f64|
                        {
                            let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                            let force = model.modified_canonical.asymptotic.weak_potential.force(&potential_distance, &potential_stiffness);
                            (model.isotensional.gibbs_free_energy_per_link(&force, &temperature) - model.isotensional.gibbs_free_energy_per_link(&force_ref, &temperature) - model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) + model.modified_canonical.asymptotic.weak_potential.gibbs_free_energy_per_link(&potential_distance_ref, &potential_stiffness, &temperature)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_potential_distance: f64|
                        {
                            let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                            let force = model.modified_canonical.asymptotic.weak_potential.force(&potential_distance, &potential_stiffness);
                            (model.isotensional.gibbs_free_energy_per_link(&force, &temperature) - model.isotensional.gibbs_free_energy_per_link(&force_ref, &temperature)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                    let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                    assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                    let potential_distance_ref = parameters.nondimensional_potential_distance_large_1*(number_of_links as f64)*link_length;
                    let residual_rel = |nondimensional_potential_stiffness|
                    {
                        let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                        let force_ref = model.modified_canonical.asymptotic.weak_potential.force(&potential_distance_ref, &potential_stiffness);
                        let integrand_numerator = |nondimensional_potential_distance: f64|
                        {
                            let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                            let force = model.modified_canonical.asymptotic.weak_potential.force(&potential_distance, &potential_stiffness);
                            (model.isotensional.relative_gibbs_free_energy(&force, &temperature) - model.isotensional.relative_gibbs_free_energy(&force_ref, &temperature) - model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature) + model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy(&potential_distance_ref, &potential_stiffness, &temperature)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_potential_distance: f64|
                        {
                            let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                            let force = model.modified_canonical.asymptotic.weak_potential.force(&potential_distance, &potential_stiffness);
                            (model.isotensional.relative_gibbs_free_energy(&force, &temperature) - model.isotensional.relative_gibbs_free_energy(&force_ref, &temperature)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                    let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                    assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                    let potential_distance_ref = parameters.nondimensional_potential_distance_large_1*(number_of_links as f64)*link_length;
                    let residual_rel = |nondimensional_potential_stiffness|
                    {
                        let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
                        let force_ref = model.modified_canonical.asymptotic.weak_potential.force(&potential_distance_ref, &potential_stiffness);
                        let integrand_numerator = |nondimensional_potential_distance: f64|
                        {
                            let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                            let force = model.modified_canonical.asymptotic.weak_potential.force(&potential_distance, &potential_stiffness);
                            (model.isotensional.relative_gibbs_free_energy_per_link(&force, &temperature) - model.isotensional.relative_gibbs_free_energy_per_link(&force_ref, &temperature) - model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature) + model.modified_canonical.asymptotic.weak_potential.relative_gibbs_free_energy_per_link(&potential_distance_ref, &potential_stiffness, &temperature)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_potential_distance: f64|
                        {
                            let potential_distance = (number_of_links as f64)*link_length*nondimensional_potential_distance;
                            let force = model.modified_canonical.asymptotic.weak_potential.force(&potential_distance, &potential_stiffness);
                            (model.isotensional.relative_gibbs_free_energy_per_link(&force, &temperature) - model.isotensional.relative_gibbs_free_energy_per_link(&force_ref, &temperature)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                    let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                    assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                    let nondimensional_potential_distance_ref = parameters.nondimensional_potential_distance_large_1;
                    let residual_rel = |nondimensional_potential_stiffness|
                    {
                        let nondimensional_force_ref = model.modified_canonical.asymptotic.weak_potential.nondimensional_force(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness);
                        let integrand_numerator = |nondimensional_potential_distance: f64|
                        {
                            let nondimensional_force = model.modified_canonical.asymptotic.weak_potential.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                            (model.isotensional.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature) - model.isotensional.nondimensional_gibbs_free_energy(&nondimensional_force_ref, &temperature) - model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) + model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness, &temperature)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_potential_distance: f64|
                        {
                            let nondimensional_force = model.modified_canonical.asymptotic.weak_potential.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                            (model.isotensional.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature) - model.isotensional.nondimensional_gibbs_free_energy(&nondimensional_force_ref, &temperature)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                    let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                    assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                    let nondimensional_potential_distance_ref = parameters.nondimensional_potential_distance_large_1;
                    let residual_rel = |nondimensional_potential_stiffness|
                    {
                        let nondimensional_force_ref = model.modified_canonical.asymptotic.weak_potential.nondimensional_force(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness);
                        let integrand_numerator = |nondimensional_potential_distance: f64|
                        {
                            let nondimensional_force = model.modified_canonical.asymptotic.weak_potential.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                            (model.isotensional.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature) - model.isotensional.nondimensional_gibbs_free_energy_per_link(&nondimensional_force_ref, &temperature) - model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature) + model.modified_canonical.asymptotic.weak_potential.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness, &temperature)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_potential_distance: f64|
                        {
                            let nondimensional_force = model.modified_canonical.asymptotic.weak_potential.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                            (model.isotensional.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature) - model.isotensional.nondimensional_gibbs_free_energy_per_link(&nondimensional_force_ref, &temperature)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                    let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                    assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                    let nondimensional_potential_distance_ref = parameters.nondimensional_potential_distance_large_1;
                    let residual_rel = |nondimensional_potential_stiffness|
                    {
                        let nondimensional_force_ref = model.modified_canonical.asymptotic.weak_potential.nondimensional_force(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness);
                        let integrand_numerator = |nondimensional_potential_distance: f64|
                        {
                            let nondimensional_force = model.modified_canonical.asymptotic.weak_potential.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                            (model.isotensional.nondimensional_relative_gibbs_free_energy(&nondimensional_force,) - model.isotensional.nondimensional_relative_gibbs_free_energy(&nondimensional_force_ref) - model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness) + model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_potential_distance: f64|
                        {
                            let nondimensional_force = model.modified_canonical.asymptotic.weak_potential.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                            (model.isotensional.nondimensional_relative_gibbs_free_energy(&nondimensional_force) - model.isotensional.nondimensional_relative_gibbs_free_energy(&nondimensional_force_ref)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                    let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                    assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
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
                    let nondimensional_potential_distance_ref = parameters.nondimensional_potential_distance_large_1;
                    let residual_rel = |nondimensional_potential_stiffness|
                    {
                        let nondimensional_force_ref = model.modified_canonical.asymptotic.weak_potential.nondimensional_force(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness);
                        let integrand_numerator = |nondimensional_potential_distance: f64|
                        {
                            let nondimensional_force = model.modified_canonical.asymptotic.weak_potential.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                            (model.isotensional.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force,) - model.isotensional.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force_ref) - model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness) + model.modified_canonical.asymptotic.weak_potential.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance_ref, &nondimensional_potential_stiffness)).powi(2)
                        };
                        let integrand_denominator = |nondimensional_potential_distance: f64|
                        {
                            let nondimensional_force = model.modified_canonical.asymptotic.weak_potential.nondimensional_force(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
                            (model.isotensional.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force) - model.isotensional.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force_ref)).powi(2)
                        };
                        let numerator = integrate(integrand_numerator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        let denominator = integrate(integrand_denominator, &parameters.nondimensional_potential_distance_large_1, &parameters.nondimensional_potential_distance_large_2, &POINTS);
                        (numerator/denominator).sqrt()
                    };
                    let residual_rel_1 = residual_rel(parameters.nondimensional_potential_stiffness_small);
                    let residual_rel_2 = residual_rel(parameters.nondimensional_potential_stiffness_small*parameters.log_log_scale);
                    let log_log_slope = -(residual_rel_2/residual_rel_1).ln()/(parameters.log_log_scale).ln();
                    assert!(residual_rel_1.abs() <= parameters.nondimensional_potential_stiffness_small);
                    assert!(residual_rel_2.abs() <= parameters.nondimensional_potential_stiffness_small/parameters.log_log_scale);
                    assert!((log_log_slope + 1.0).abs() <= parameters.log_log_tol);
                }
            }
        }
    }
}
