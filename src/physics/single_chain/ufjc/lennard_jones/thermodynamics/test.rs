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
        let _ = LENNARDJONESFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference);
    }
    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            assert_eq!(number_of_links, LENNARDJONESFJC::init(number_of_links, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).number_of_links);
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
            assert_eq!(link_length, LENNARDJONESFJC::init(parameters.number_of_links_minimum, link_length, parameters.hinge_mass_reference, parameters.link_stiffness_reference).link_length);
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
            assert_eq!(hinge_mass, LENNARDJONESFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, hinge_mass, parameters.link_stiffness_reference).hinge_mass);
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
            assert_eq!(link_stiffness, LENNARDJONESFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, link_stiffness).link_stiffness);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            assert_eq!(number_of_links, model.number_of_links);
            assert_eq!(link_length, model.link_length);
            assert_eq!(hinge_mass, model.hinge_mass);
            assert_eq!(link_stiffness, model.link_stiffness);
        }
    }
}
mod legendre_asymptotic
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.asymptotic.end_to_end_length(&force, &temperature);
            let force_out = model.isometric.asymptotic.legendre.force(&end_to_end_length, &temperature);
            let residual_abs = &force - &force_out;
            let residual_rel = &residual_abs/&force;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_end_to_end_length_per_link= model.isotensional.asymptotic.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature);
            let nondimensional_force_out = model.isometric.asymptotic.legendre.nondimensional_force(&nondimensional_end_to_end_length_per_link, &temperature);
            let residual_abs = &nondimensional_force - &nondimensional_force_out;
            let residual_rel = &residual_abs/&nondimensional_force;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.asymptotic.end_to_end_length(&force, &temperature);
            let helmholtz_free_energy_legendre = model.isotensional.asymptotic.gibbs_free_energy(&force, &temperature) + force*end_to_end_length;
            let helmholtz_free_energy_legendre_out = model.isometric.asymptotic.legendre.helmholtz_free_energy(&end_to_end_length, &temperature);
            let residual_abs = &helmholtz_free_energy_legendre - &helmholtz_free_energy_legendre_out + BOLTZMANN_CONSTANT*temperature*(0.5*(2.0*PI*BOLTZMANN_CONSTANT*temperature/link_stiffness).ln() + (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln());
            let residual_rel = &residual_abs/&helmholtz_free_energy_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.asymptotic.end_to_end_length(&force, &temperature);
            let end_to_end_length_per_link = model.isotensional.asymptotic.end_to_end_length_per_link(&force, &temperature);
            let helmholtz_free_energy_per_link_legendre = model.isotensional.asymptotic.gibbs_free_energy_per_link(&force, &temperature) + force*end_to_end_length_per_link;
            let helmholtz_free_energy_per_link_legendre_out = model.isometric.asymptotic.legendre.helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
            let residual_abs = &helmholtz_free_energy_per_link_legendre - &helmholtz_free_energy_per_link_legendre_out + BOLTZMANN_CONSTANT*temperature*(0.5*(2.0*PI*BOLTZMANN_CONSTANT*temperature/link_stiffness).ln() + (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())/(number_of_links as f64);
            let residual_rel = &residual_abs/&helmholtz_free_energy_per_link_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.asymptotic.end_to_end_length(&force, &temperature);
            let relative_helmholtz_free_energy_legendre = model.isotensional.asymptotic.relative_gibbs_free_energy(&force, &temperature) + force*end_to_end_length;
            let relative_helmholtz_free_energy_legendre_out = model.isometric.asymptotic.legendre.relative_helmholtz_free_energy(&end_to_end_length, &temperature);
            let residual_abs = &relative_helmholtz_free_energy_legendre - &relative_helmholtz_free_energy_legendre_out;
            let residual_rel = &residual_abs/&relative_helmholtz_free_energy_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.asymptotic.end_to_end_length(&force, &temperature);
            let end_to_end_length_per_link = model.isotensional.asymptotic.end_to_end_length_per_link(&force, &temperature);
            let relative_helmholtz_free_energy_per_link_legendre = model.isotensional.asymptotic.relative_gibbs_free_energy_per_link(&force, &temperature) + force*end_to_end_length_per_link;
            let relative_helmholtz_free_energy_per_link_legendre_out = model.isometric.asymptotic.legendre.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
            let residual_abs = &relative_helmholtz_free_energy_per_link_legendre - &relative_helmholtz_free_energy_per_link_legendre_out;
            let residual_rel = &residual_abs/&relative_helmholtz_free_energy_per_link_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_end_to_end_length = model.isotensional.asymptotic.nondimensional_end_to_end_length(&nondimensional_force, &temperature);
            let nondimensional_end_to_end_length_per_link = model.isotensional.asymptotic.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature);
            let nondimensional_helmholtz_free_energy_legendre = model.isotensional.asymptotic.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature) + nondimensional_force*nondimensional_end_to_end_length;
            let nondimensional_helmholtz_free_energy_legendre_out = model.isometric.asymptotic.legendre.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let residual_abs = &nondimensional_helmholtz_free_energy_legendre - &nondimensional_helmholtz_free_energy_legendre_out + 0.5*(2.0*PI*BOLTZMANN_CONSTANT*temperature/link_stiffness).ln() + (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln();
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_end_to_end_length_per_link = model.isotensional.asymptotic.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature);
            let nondimensional_helmholtz_free_energy_per_link_legendre = model.isotensional.asymptotic.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature) + nondimensional_force*nondimensional_end_to_end_length_per_link;
            let nondimensional_helmholtz_free_energy_per_link_legendre_out = model.isometric.asymptotic.legendre.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
            let residual_abs = &nondimensional_helmholtz_free_energy_per_link_legendre - &nondimensional_helmholtz_free_energy_per_link_legendre_out + (0.5*(2.0*PI*BOLTZMANN_CONSTANT*temperature/link_stiffness).ln() + (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())/(number_of_links as f64);
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_end_to_end_length = model.isotensional.asymptotic.nondimensional_end_to_end_length(&nondimensional_force, &temperature);
            let nondimensional_end_to_end_length_per_link = model.isotensional.asymptotic.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature);
            let nondimensional_relative_helmholtz_free_energy_legendre = model.isotensional.asymptotic.nondimensional_relative_gibbs_free_energy(&nondimensional_force, &temperature) + nondimensional_force*nondimensional_end_to_end_length;
            let nondimensional_relative_helmholtz_free_energy_legendre_out = model.isometric.asymptotic.legendre.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let residual_abs = &nondimensional_relative_helmholtz_free_energy_legendre - &nondimensional_relative_helmholtz_free_energy_legendre_out;
            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_end_to_end_length_per_link = model.isotensional.asymptotic.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature);
            let nondimensional_relative_helmholtz_free_energy_per_link_legendre = model.isotensional.asymptotic.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force, &temperature) + nondimensional_force*nondimensional_end_to_end_length_per_link;
            let nondimensional_relative_helmholtz_free_energy_per_link_legendre_out = model.isometric.asymptotic.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
            let residual_abs = &nondimensional_relative_helmholtz_free_energy_per_link_legendre - &nondimensional_relative_helmholtz_free_energy_per_link_legendre_out;
            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
        }
    }
}
mod legendre_asymptotic_reduced
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.asymptotic.reduced.end_to_end_length(&force, &temperature);
            let force_out = model.isometric.asymptotic.reduced.legendre.force(&end_to_end_length, &temperature);
            let residual_abs = &force - &force_out;
            let residual_rel = &residual_abs/&force;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_end_to_end_length_per_link= model.isotensional.asymptotic.reduced.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature);
            let nondimensional_force_out = model.isometric.asymptotic.reduced.legendre.nondimensional_force(&nondimensional_end_to_end_length_per_link, &temperature);
            let residual_abs = &nondimensional_force - &nondimensional_force_out;
            let residual_rel = &residual_abs/&nondimensional_force;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.asymptotic.reduced.end_to_end_length(&force, &temperature);
            let helmholtz_free_energy_legendre = model.isotensional.asymptotic.reduced.gibbs_free_energy(&force, &temperature) + force*end_to_end_length;
            let helmholtz_free_energy_legendre_out = model.isometric.asymptotic.reduced.legendre.helmholtz_free_energy(&end_to_end_length, &temperature);
            let residual_abs = &helmholtz_free_energy_legendre - &helmholtz_free_energy_legendre_out + BOLTZMANN_CONSTANT*temperature*(0.5*(2.0*PI*BOLTZMANN_CONSTANT*temperature/link_stiffness).ln() + (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln());
            let residual_rel = &residual_abs/&helmholtz_free_energy_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.asymptotic.reduced.end_to_end_length(&force, &temperature);
            let end_to_end_length_per_link = model.isotensional.asymptotic.reduced.end_to_end_length_per_link(&force, &temperature);
            let helmholtz_free_energy_per_link_legendre = model.isotensional.asymptotic.reduced.gibbs_free_energy_per_link(&force, &temperature) + force*end_to_end_length_per_link;
            let helmholtz_free_energy_per_link_legendre_out = model.isometric.asymptotic.reduced.legendre.helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
            let residual_abs = &helmholtz_free_energy_per_link_legendre - &helmholtz_free_energy_per_link_legendre_out + BOLTZMANN_CONSTANT*temperature*(0.5*(2.0*PI*BOLTZMANN_CONSTANT*temperature/link_stiffness).ln() + (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())/(number_of_links as f64);
            let residual_rel = &residual_abs/&helmholtz_free_energy_per_link_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.asymptotic.reduced.end_to_end_length(&force, &temperature);
            let relative_helmholtz_free_energy_legendre = model.isotensional.asymptotic.reduced.relative_gibbs_free_energy(&force, &temperature) + force*end_to_end_length;
            let relative_helmholtz_free_energy_legendre_out = model.isometric.asymptotic.reduced.legendre.relative_helmholtz_free_energy(&end_to_end_length, &temperature);
            let residual_abs = &relative_helmholtz_free_energy_legendre - &relative_helmholtz_free_energy_legendre_out;
            let residual_rel = &residual_abs/&relative_helmholtz_free_energy_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let end_to_end_length = model.isotensional.asymptotic.reduced.end_to_end_length(&force, &temperature);
            let end_to_end_length_per_link = model.isotensional.asymptotic.reduced.end_to_end_length_per_link(&force, &temperature);
            let relative_helmholtz_free_energy_per_link_legendre = model.isotensional.asymptotic.reduced.relative_gibbs_free_energy_per_link(&force, &temperature) + force*end_to_end_length_per_link;
            let relative_helmholtz_free_energy_per_link_legendre_out = model.isometric.asymptotic.reduced.legendre.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
            let residual_abs = &relative_helmholtz_free_energy_per_link_legendre - &relative_helmholtz_free_energy_per_link_legendre_out;
            let residual_rel = &residual_abs/&relative_helmholtz_free_energy_per_link_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_end_to_end_length = model.isotensional.asymptotic.reduced.nondimensional_end_to_end_length(&nondimensional_force, &temperature);
            let nondimensional_end_to_end_length_per_link = model.isotensional.asymptotic.reduced.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature);
            let nondimensional_helmholtz_free_energy_legendre = model.isotensional.asymptotic.reduced.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature) + nondimensional_force*nondimensional_end_to_end_length;
            let nondimensional_helmholtz_free_energy_legendre_out = model.isometric.asymptotic.reduced.legendre.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let residual_abs = &nondimensional_helmholtz_free_energy_legendre - &nondimensional_helmholtz_free_energy_legendre_out + 0.5*(2.0*PI*BOLTZMANN_CONSTANT*temperature/link_stiffness).ln() + (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln();
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_end_to_end_length_per_link = model.isotensional.asymptotic.reduced.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature);
            let nondimensional_helmholtz_free_energy_per_link_legendre = model.isotensional.asymptotic.reduced.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature) + nondimensional_force*nondimensional_end_to_end_length_per_link;
            let nondimensional_helmholtz_free_energy_per_link_legendre_out = model.isometric.asymptotic.reduced.legendre.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
            let residual_abs = &nondimensional_helmholtz_free_energy_per_link_legendre - &nondimensional_helmholtz_free_energy_per_link_legendre_out + (0.5*(2.0*PI*BOLTZMANN_CONSTANT*temperature/link_stiffness).ln() + (8.0*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln())/(number_of_links as f64);
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_end_to_end_length = model.isotensional.asymptotic.reduced.nondimensional_end_to_end_length(&nondimensional_force, &temperature);
            let nondimensional_end_to_end_length_per_link = model.isotensional.asymptotic.reduced.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature);
            let nondimensional_relative_helmholtz_free_energy_legendre = model.isotensional.asymptotic.reduced.nondimensional_relative_gibbs_free_energy(&nondimensional_force, &temperature) + nondimensional_force*nondimensional_end_to_end_length;
            let nondimensional_relative_helmholtz_free_energy_legendre_out = model.isometric.asymptotic.reduced.legendre.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let residual_abs = &nondimensional_relative_helmholtz_free_energy_legendre - &nondimensional_relative_helmholtz_free_energy_legendre_out;
            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
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
            let model = LENNARDJONESFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let lambda_max = (13.0/7.0 as f64).powf(1.0/6.0);
            let nondimensional_force_max = link_stiffness/BOLTZMANN_CONSTANT/temperature*link_length.powi(2)/6.0*(lambda_max.powi(-7) - lambda_max.powi(-13));
            let nondimensional_force = nondimensional_force_max*rng.gen::<f64>();
            let nondimensional_end_to_end_length_per_link = model.isotensional.asymptotic.reduced.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature);
            let nondimensional_relative_helmholtz_free_energy_per_link_legendre = model.isotensional.asymptotic.reduced.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force, &temperature) + nondimensional_force*nondimensional_end_to_end_length_per_link;
            let nondimensional_relative_helmholtz_free_energy_per_link_legendre_out = model.isometric.asymptotic.reduced.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
            let residual_abs = &nondimensional_relative_helmholtz_free_energy_per_link_legendre - &nondimensional_relative_helmholtz_free_energy_per_link_legendre_out;
            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link_legendre;
            assert!(residual_abs.abs() <= parameters.abs_tol || residual_rel.abs() <= parameters.rel_tol);
        }
    }
}