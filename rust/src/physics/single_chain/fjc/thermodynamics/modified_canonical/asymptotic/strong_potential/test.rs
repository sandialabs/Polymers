#![cfg(test)]
use super::*;
use crate::physics::single_chain::test::Parameters as DefaultParameters;
pub struct Parameters
{
    pub abs_tol: f64,
    pub rel_tol: f64,
    pub number_of_loops: u32,
    pub hinge_mass_reference: f64,
    pub hinge_mass_scale: f64,
    pub link_length_reference: f64,
    pub link_length_scale: f64,
    pub number_of_links_minimum: u8,
    pub number_of_links_maximum: u8,
    pub nondimensional_potential_distance_small: f64,
    pub nondimensional_potential_stiffness_reference: f64,
    pub nondimensional_potential_stiffness_scale: f64,
    pub temperature_reference: f64,
    pub temperature_scale: f64,
}
impl Default for Parameters
{
    fn default() -> Self
    {
        Self
        {
            number_of_loops: 888,
            abs_tol: DefaultParameters::default().abs_tol,
            rel_tol: DefaultParameters::default().rel_tol,
            hinge_mass_reference: DefaultParameters::default().hinge_mass_reference,
            hinge_mass_scale: DefaultParameters::default().hinge_mass_scale,
            link_length_reference: DefaultParameters::default().link_length_reference,
            link_length_scale: DefaultParameters::default().link_length_scale,
            number_of_links_minimum: DefaultParameters::default().number_of_links_minimum,
            number_of_links_maximum: DefaultParameters::default().number_of_links_maximum,
            nondimensional_potential_distance_small: DefaultParameters::default().nondimensional_potential_distance_small,
            nondimensional_potential_stiffness_reference: DefaultParameters::default().nondimensional_potential_stiffness_reference,
            nondimensional_potential_stiffness_scale: DefaultParameters::default().nondimensional_potential_stiffness_scale,
            temperature_reference: DefaultParameters::default().temperature_reference,
            temperature_scale: DefaultParameters::default().temperature_scale,
        }
    }
}
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
    #[ignore]
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
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let force = model.force(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = &force/BOLTZMANN_CONSTANT/temperature*link_length - &nondimensional_force;
            let residual_rel = &residual_abs/&nondimensional_force;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    #[ignore]
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
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let helmholtz_free_energy = model.helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = helmholtz_free_energy/BOLTZMANN_CONSTANT/temperature - nondimensional_helmholtz_free_energy;
            let residual_rel = residual_abs/nondimensional_helmholtz_free_energy;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    #[ignore]
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
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let residual_abs = helmholtz_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - nondimensional_helmholtz_free_energy_per_link;
            let residual_rel = residual_abs/nondimensional_helmholtz_free_energy_per_link;
            assert!(residual_abs.abs() <= parameters.abs_tol);
            assert!(residual_rel.abs() <= parameters.rel_tol);
        }
    }
    #[test]
    #[ignore]
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
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_small*(0.5 + (0.5 - rng.gen::<f64>()));
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
    #[ignore]
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
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_small*(0.5 + (0.5 - rng.gen::<f64>()));
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
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_small*(0.5 + (0.5 - rng.gen::<f64>()));
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
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_small*(0.5 + (0.5 - rng.gen::<f64>()));
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
            let residual_abs = nondimensional_relative_helmholtz_free_energy/(number_of_links as f64) - nondimensional_relative_helmholtz_free_energy_per_link;
            let residual_rel = residual_abs/nondimensional_relative_helmholtz_free_energy_per_link;
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
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_small*(0.5 + (0.5 - rng.gen::<f64>()));
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let helmholtz_free_energy = model.helmholtz_free_energy(&potential_distance, &potential_stiffness, &temperature);
            let helmholtz_free_energy_0 = model.helmholtz_free_energy(&(ZERO*(number_of_links as f64)*link_length), &potential_stiffness, &temperature);
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
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_small*(0.5 + (0.5 - rng.gen::<f64>()));
            let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
            let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powi(2)*BOLTZMANN_CONSTANT*temperature;
            let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
            let helmholtz_free_energy_per_link_0 = model.helmholtz_free_energy_per_link(&(ZERO*(number_of_links as f64)*link_length), &potential_stiffness, &temperature);
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
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_small*(0.5 + (0.5 - rng.gen::<f64>()));
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
            let nondimensional_potential_distance = parameters.nondimensional_potential_distance_small*(0.5 + (0.5 - rng.gen::<f64>()));
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
}