#![cfg(test)]
use super::*;
use crate::physics::
{
    BOLTZMANN_CONSTANT,
    PLANCK_CONSTANT
};
use crate::physics::single_chain::test::Parameters;
use std::f64::consts::PI;
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let force = model.isometric.force(&end_to_end_length, &temperature);
            let helmholtz_free_energy = model.isometric.helmholtz_free_energy(&end_to_end_length, &temperature);
            let helmholtz_free_energy_out = model.isotensional.gibbs_free_energy(&force, &temperature) + force*end_to_end_length;
            let residual_abs = &helmholtz_free_energy - &helmholtz_free_energy_out - BOLTZMANN_CONSTANT*temperature*(4.0*(-(number_of_links as f64)*link_length/persistance_length).exp().acos().sin()*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln();
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let end_to_end_length_per_link = nondimensional_end_to_end_length_per_link*link_length;
            let force = model.isometric.force(&end_to_end_length, &temperature);
            let helmholtz_free_energy_per_link = model.isometric.helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
            let helmholtz_free_energy_per_link_out = model.isotensional.gibbs_free_energy_per_link(&force, &temperature) + force*end_to_end_length_per_link;
            let residual_abs = &helmholtz_free_energy_per_link - &helmholtz_free_energy_per_link_out - BOLTZMANN_CONSTANT*temperature*(4.0*(-(number_of_links as f64)*link_length/persistance_length).exp().acos().sin()*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()/(number_of_links as f64);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64);
            let nondimensional_force = model.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link);
            let nondimensional_helmholtz_free_energy = model.isometric.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let nondimensional_helmholtz_free_energy_out = model.isotensional.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature) + nondimensional_force*nondimensional_end_to_end_length;
            let residual_abs = &nondimensional_helmholtz_free_energy - &nondimensional_helmholtz_free_energy_out - (4.0*(-(number_of_links as f64)*link_length/persistance_length).exp().acos().sin()*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln();
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force = model.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link);
            let nondimensional_helmholtz_free_energy_per_link = model.isometric.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
            let nondimensional_helmholtz_free_energy_per_link_out = model.isotensional.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature) + nondimensional_force*nondimensional_end_to_end_length_per_link;
            let residual_abs = &nondimensional_helmholtz_free_energy_per_link - &nondimensional_helmholtz_free_energy_per_link_out - (4.0*(-(number_of_links as f64)*link_length/persistance_length).exp().acos().sin()*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()/(number_of_links as f64);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.end_to_end_length(&force, &temperature);
            let gibbs_free_energy = model.isotensional.gibbs_free_energy(&force, &temperature);
            let gibbs_free_energy_out = model.isometric.helmholtz_free_energy(&end_to_end_length, &temperature) - force*end_to_end_length;
            let residual_abs = &gibbs_free_energy - &gibbs_free_energy_out + BOLTZMANN_CONSTANT*temperature*(4.0*(-(number_of_links as f64)*link_length/persistance_length).exp().acos().sin()*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln();
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.end_to_end_length(&force, &temperature);
            let end_to_end_length_per_link = end_to_end_length/(number_of_links as f64);
            let gibbs_free_energy_per_link = model.isotensional.gibbs_free_energy_per_link(&force, &temperature);
            let gibbs_free_energy_per_link_out = model.isometric.helmholtz_free_energy_per_link(&end_to_end_length, &temperature) - force*end_to_end_length_per_link;
            let residual_abs = &gibbs_free_energy_per_link - &gibbs_free_energy_per_link_out + BOLTZMANN_CONSTANT*temperature*(4.0*(-(number_of_links as f64)*link_length/persistance_length).exp().acos().sin()*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()/(number_of_links as f64);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = model.isotensional.nondimensional_end_to_end_length(&nondimensional_force);
            let nondimensional_end_to_end_length_per_link = nondimensional_end_to_end_length/(number_of_links as f64);
            let nondimensional_gibbs_free_energy = model.isotensional.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature);
            let nondimensional_gibbs_free_energy_out = model.isometric.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature) - nondimensional_force*nondimensional_end_to_end_length;
            let residual_abs = &nondimensional_gibbs_free_energy - &nondimensional_gibbs_free_energy_out + (4.0*(-(number_of_links as f64)*link_length/persistance_length).exp().acos().sin()*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln();
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length_per_link = model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
            let nondimensional_gibbs_free_energy_per_link = model.isotensional.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature);
            let nondimensional_gibbs_free_energy_per_link_out = model.isometric.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature) - nondimensional_force*nondimensional_end_to_end_length_per_link;
            let residual_abs = &nondimensional_gibbs_free_energy_per_link - &nondimensional_gibbs_free_energy_per_link_out + (4.0*(-(number_of_links as f64)*link_length/persistance_length).exp().acos().sin()*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()/(number_of_links as f64);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
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
            let number_of_links: u8 = parameters.number_of_links_minimum;
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let persistance_length = parameters.nondimensional_persistance_length_small*(number_of_links as f64)*link_length;
            let model = WLC::init(number_of_links, link_length, hinge_mass, persistance_length);
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