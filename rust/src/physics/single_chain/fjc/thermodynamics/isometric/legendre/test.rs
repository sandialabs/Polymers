#![cfg(test)]
use super::*;
use crate::physics::single_chain::fjc::thermodynamics::isometric::test::Parameters;
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
mod normalization
{
    use super::*;
    use rand::Rng;
    use crate::physics::single_chain::test::integrate;
    #[test]
    fn equilibrium_distribution()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let integrand = |end_to_end_length: f64| 4.0*PI*end_to_end_length.powi(2)*model.equilibrium_distribution(&end_to_end_length);
            let integral = integrate(integrand, &ZERO, &(ONE*model.contour_length), &POINTS);
            assert!((integral - 1.0).abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn nondimensional_equilibrium_distribution()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let integrand = |nondimensional_end_to_end_length_per_link_per_link: f64| 4.0*PI*nondimensional_end_to_end_length_per_link_per_link.powi(2)*model.nondimensional_equilibrium_distribution(&nondimensional_end_to_end_length_per_link_per_link);
            let integral = integrate(integrand, &ZERO, &ONE, &POINTS);
            assert!((integral - 1.0).abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn equilibrium_radial_distribution()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let integrand = |end_to_end_length: f64| model.equilibrium_radial_distribution(&end_to_end_length);
            let integral = integrate(integrand, &ZERO, &(ONE*model.contour_length), &POINTS);
            assert!((integral - 1.0).abs() <= parameters.rel_tol);
        }
    }
    #[test]
    fn nondimensional_equilibrium_radial_distribution()
    {
        let parameters = Parameters::default();
        let mut rng = rand::thread_rng();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let integrand = |nondimensional_end_to_end_length_per_link_per_link: f64| model.nondimensional_equilibrium_radial_distribution(&nondimensional_end_to_end_length_per_link_per_link);
            let integral = integrate(integrand, &ZERO, &ONE, &POINTS);
            assert!((integral - 1.0).abs() <= parameters.rel_tol);
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
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_force = model.nondimensional_force(&nondimensional_end_to_end_length_per_link);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let helmholtz_free_energy = model.helmholtz_free_energy(&end_to_end_length, &temperature);
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
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
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
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link);
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&end_to_end_length, &temperature);
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
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
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
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let helmholtz_free_energy = model.helmholtz_free_energy(&end_to_end_length, &temperature);
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
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
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&end_to_end_length, &temperature);
            let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
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
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
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
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link);
            let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
            let residual_abs = nondimensional_relative_helmholtz_free_energy/(number_of_links as f64) - nondimensional_relative_helmholtz_free_energy_per_link;
            let residual_rel = residual_abs/nondimensional_relative_helmholtz_free_energy_per_link;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
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
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let helmholtz_free_energy = model.helmholtz_free_energy(&end_to_end_length, &temperature);
            let helmholtz_free_energy_0 = model.helmholtz_free_energy(&ZERO, &temperature);
            let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&end_to_end_length, &temperature);
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
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
            let helmholtz_free_energy_per_link_0 = model.helmholtz_free_energy_per_link(&ZERO, &temperature);
            let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
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
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
            let nondimensional_helmholtz_free_energy_0 = model.nondimensional_helmholtz_free_energy(&ZERO, &temperature);
            let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link);
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
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
            let nondimensional_helmholtz_free_energy_per_link_0 = model.nondimensional_helmholtz_free_energy_per_link(&ZERO, &temperature);
            let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
            let residual_abs = &nondimensional_helmholtz_free_energy_per_link - &nondimensional_helmholtz_free_energy_per_link_0 - &nondimensional_relative_helmholtz_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link_0;
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let gibbs_free_energy = model.gibbs_free_energy(&end_to_end_length, &temperature);
            let gibbs_free_energy_0 = model.gibbs_free_energy(&ZERO, &temperature);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
            let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&end_to_end_length, &temperature);
            let gibbs_free_energy_per_link_0 = model.gibbs_free_energy_per_link(&ZERO, &temperature);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
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
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force_0 = model.force(&(ZERO*(number_of_links as f64)*link_length), &temperature);
            assert!(force_0.abs() <= 3.0*BOLTZMANN_CONSTANT*temperature/link_length*ZERO);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_force_0 = model.nondimensional_force(&ZERO);
            assert!(nondimensional_force_0.abs() <= 3.0*ZERO);
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
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_helmholtz_free_energy_0 = model.relative_helmholtz_free_energy(&(ZERO*(number_of_links as f64)*link_length), &temperature);
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
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_helmholtz_free_energy_per_link_0 = model.relative_helmholtz_free_energy_per_link(&(ZERO*(number_of_links as f64)*link_length), &temperature);
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
            let nondimensional_relative_helmholtz_free_energy_0 = model.nondimensional_relative_helmholtz_free_energy(&ZERO);
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
            let nondimensional_relative_helmholtz_free_energy_per_link_0 = model.nondimensional_relative_helmholtz_free_energy_per_link(&ZERO);
            assert!(nondimensional_relative_helmholtz_free_energy_per_link_0.abs() <= ZERO);
        }
    }
    #[test]
    fn equilibrium_radial_distribution()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let equilibrium_radial_distribution_0 = model.equilibrium_radial_distribution(&(ZERO*(number_of_links as f64)*link_length));
            assert!(equilibrium_radial_distribution_0.abs() <= ZERO);
        }
    }
    #[test]
    fn nondimensional_equilibrium_radial_distribution()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_equilibrium_radial_distribution_0 = model.nondimensional_equilibrium_radial_distribution(&ZERO);
            assert!(nondimensional_equilibrium_radial_distribution_0.abs() <= ZERO);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
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
            let model = FJC::init(number_of_links, link_length, hinge_mass);
            let nondimensional_relative_gibbs_free_energy_per_link_0 = model.nondimensional_relative_gibbs_free_energy_per_link(&ZERO);
            assert!(nondimensional_relative_gibbs_free_energy_per_link_0.abs() <= ZERO);
        }
    }
}
