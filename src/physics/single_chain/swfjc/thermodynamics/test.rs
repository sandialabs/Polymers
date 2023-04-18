#![cfg(test)]
use super::*;
use std::f64::consts::PI;
use crate::physics::
{
    BOLTZMANN_CONSTANT,
    PLANCK_CONSTANT
};
use crate::physics::single_chain::test::Parameters;
mod base
{
    use super::*;
    use rand::Rng;
    #[test]
    fn init()
    {
        let parameters = Parameters::default();
        let _ = SWFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.well_width_reference);
    }
    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            assert_eq!(number_of_links, SWFJC::init(number_of_links, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.well_width_reference).number_of_links);
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
            assert_eq!(link_length, SWFJC::init(parameters.number_of_links_minimum, link_length, parameters.hinge_mass_reference, parameters.well_width_reference).link_length);
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
            assert_eq!(hinge_mass, SWFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, hinge_mass, parameters.well_width_reference).hinge_mass);
        }
    }
    #[test]
    fn well_width()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(well_width, SWFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, well_width).well_width);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            assert_eq!(number_of_links, model.number_of_links);
            assert_eq!(link_length, model.link_length);
            assert_eq!(hinge_mass, model.hinge_mass);
            assert_eq!(well_width, model.well_width);
        }
    }
}
mod legendre
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.end_to_end_length(&force, &temperature);
            let force_out = model.isometric.legendre.force(&end_to_end_length, &temperature);
            let residual_abs = &force - &force_out;
            let residual_rel = &residual_abs/&force;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length_per_link= model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
            let nondimensional_force_out = model.isometric.legendre.nondimensional_force(&nondimensional_end_to_end_length_per_link);
            let residual_abs = &nondimensional_force - &nondimensional_force_out;
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.end_to_end_length(&force, &temperature);
            let helmholtz_free_energy_legendre = model.isotensional.gibbs_free_energy(&force, &temperature) + force*end_to_end_length;
            let helmholtz_free_energy_legendre_out = model.isometric.legendre.helmholtz_free_energy(&end_to_end_length, &temperature);
            let residual_abs = &helmholtz_free_energy_legendre - &helmholtz_free_energy_legendre_out + BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln();
            let residual_rel = &residual_abs/&helmholtz_free_energy_legendre;
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.end_to_end_length(&force, &temperature);
            let end_to_end_length_per_link = model.isotensional.end_to_end_length_per_link(&force, &temperature);
            let helmholtz_free_energy_per_link_legendre = model.isotensional.gibbs_free_energy_per_link(&force, &temperature) + force*end_to_end_length_per_link;
            let helmholtz_free_energy_per_link_legendre_out = model.isometric.legendre.helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
            let residual_abs = &helmholtz_free_energy_per_link_legendre - &helmholtz_free_energy_per_link_legendre_out + BOLTZMANN_CONSTANT*temperature*(8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()/(number_of_links as f64);
            let residual_rel = &residual_abs/&helmholtz_free_energy_per_link_legendre;
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.end_to_end_length(&force, &temperature);
            let relative_helmholtz_free_energy_legendre = model.isotensional.relative_gibbs_free_energy(&force, &temperature) + force*end_to_end_length;
            let relative_helmholtz_free_energy_legendre_out = model.isometric.legendre.relative_helmholtz_free_energy(&end_to_end_length, &temperature);
            let residual_abs = &relative_helmholtz_free_energy_legendre - &relative_helmholtz_free_energy_legendre_out;
            let residual_rel = &residual_abs/&relative_helmholtz_free_energy_legendre;
            assert!(residual_rel.abs() <= 3e1 * parameters.rel_tol);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.end_to_end_length(&force, &temperature);
            let end_to_end_length_per_link = model.isotensional.end_to_end_length_per_link(&force, &temperature);
            let relative_helmholtz_free_energy_per_link_legendre = model.isotensional.relative_gibbs_free_energy_per_link(&force, &temperature) + force*end_to_end_length_per_link;
            let relative_helmholtz_free_energy_per_link_legendre_out = model.isometric.legendre.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
            let residual_abs = &relative_helmholtz_free_energy_per_link_legendre - &relative_helmholtz_free_energy_per_link_legendre_out;
            let residual_rel = &residual_abs/&relative_helmholtz_free_energy_per_link_legendre;
            assert!(residual_rel.abs() <= 3e1 * parameters.rel_tol);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = model.isotensional.nondimensional_end_to_end_length(&nondimensional_force);
            let nondimensional_end_to_end_length_per_link = model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
            let nondimensional_helmholtz_free_energy_legendre = model.isotensional.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature) + nondimensional_force*nondimensional_end_to_end_length;
            let nondimensional_helmholtz_free_energy_legendre_out = model.isometric.legendre.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let residual_abs = &nondimensional_helmholtz_free_energy_legendre - &nondimensional_helmholtz_free_energy_legendre_out + (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln();
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length_per_link = model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
            let nondimensional_helmholtz_free_energy_per_link_legendre = model.isotensional.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature) + nondimensional_force*nondimensional_end_to_end_length_per_link;
            let nondimensional_helmholtz_free_energy_per_link_legendre_out = model.isometric.legendre.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
            let residual_abs = &nondimensional_helmholtz_free_energy_per_link_legendre - &nondimensional_helmholtz_free_energy_per_link_legendre_out + (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()/(number_of_links as f64);
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link_legendre;
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length = model.isotensional.nondimensional_end_to_end_length(&nondimensional_force);
            let nondimensional_end_to_end_length_per_link = model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
            let nondimensional_relative_helmholtz_free_energy_legendre = model.isotensional.nondimensional_relative_gibbs_free_energy(&nondimensional_force) + nondimensional_force*nondimensional_end_to_end_length;
            let nondimensional_relative_helmholtz_free_energy_legendre_out = model.isometric.legendre.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link);
            let residual_abs = &nondimensional_relative_helmholtz_free_energy_legendre - &nondimensional_relative_helmholtz_free_energy_legendre_out;
            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_legendre;
            assert!(residual_rel.abs() <= 3e1 * parameters.rel_tol);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_end_to_end_length_per_link = model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
            let nondimensional_relative_helmholtz_free_energy_per_link_legendre = model.isotensional.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force) + nondimensional_force*nondimensional_end_to_end_length_per_link;
            let nondimensional_relative_helmholtz_free_energy_per_link_legendre_out = model.isometric.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
            let residual_abs = &nondimensional_relative_helmholtz_free_energy_per_link_legendre - &nondimensional_relative_helmholtz_free_energy_per_link_legendre_out;
            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link_legendre;
            assert!(residual_rel.abs() <= 3e1 * parameters.rel_tol);
        }
    }
}