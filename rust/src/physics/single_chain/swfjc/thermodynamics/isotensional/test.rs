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
    pub well_width_reference: f64,
    pub well_width_scale: f64,
    pub nondimensional_force_reference: f64,
    pub nondimensional_force_scale: f64,
    pub temperature_reference: f64,
    pub temperature_scale: f64,
}
impl Default for Parameters
{
    fn default() -> Self
    {
        Self
        {
            number_of_loops: 8888,
            abs_tol: DefaultParameters::default().abs_tol,
            rel_tol: DefaultParameters::default().rel_tol,
            hinge_mass_reference: DefaultParameters::default().hinge_mass_reference,
            hinge_mass_scale: DefaultParameters::default().hinge_mass_scale,
            link_length_reference: DefaultParameters::default().link_length_reference,
            link_length_scale: DefaultParameters::default().link_length_scale,
            number_of_links_minimum: DefaultParameters::default().number_of_links_minimum,
            number_of_links_maximum: DefaultParameters::default().number_of_links_maximum,
            well_width_reference: DefaultParameters::default().well_width_reference,
            well_width_scale: DefaultParameters::default().well_width_scale,
            nondimensional_force_reference: DefaultParameters::default().nondimensional_force_reference,
            nondimensional_force_scale: DefaultParameters::default().nondimensional_force_scale,
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
    fn number_of_links_and_link_length_and_hinge_mass_and_well_width()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_length = rng.gen::<f64>();
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(link_length, SWFJC::init(number_of_links, link_length, hinge_mass, well_width).link_length);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let gibbs_free_energy = model.gibbs_free_energy(&force, &temperature);
            let gibbs_free_energy_0 = model.gibbs_free_energy(&(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), &temperature);
            let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&force, &temperature);
            let residual_abs = &gibbs_free_energy - &gibbs_free_energy_0 - &relative_gibbs_free_energy;
            let residual_rel = &residual_abs/&gibbs_free_energy_0;
            assert!(residual_abs.abs() <= parameters.abs_tol);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
            let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&force, &temperature);
            let gibbs_free_energy_per_link_0 = model.gibbs_free_energy_per_link(&(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), &temperature);
            let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&force, &temperature);
            let residual_abs = &gibbs_free_energy_per_link - &gibbs_free_energy_per_link_0 - &relative_gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&gibbs_free_energy_per_link_0;
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature);
            let nondimensional_gibbs_free_energy_0 = model.nondimensional_gibbs_free_energy(&ZERO, &temperature);
            let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
            let residual_abs = &nondimensional_gibbs_free_energy - &nondimensional_gibbs_free_energy_0 - &nondimensional_relative_gibbs_free_energy;
            let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_0;
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature);
            let nondimensional_gibbs_free_energy_per_link_0 = model.nondimensional_gibbs_free_energy_per_link(&ZERO, &temperature);
            let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
            let residual_abs = &nondimensional_gibbs_free_energy_per_link - &nondimensional_gibbs_free_energy_per_link_0 - &nondimensional_relative_gibbs_free_energy_per_link;
            let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link_0;
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
    fn relative_gibbs_free_energy()
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
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_gibbs_free_energy_0 = model.relative_gibbs_free_energy(&(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), &temperature);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
            let relative_gibbs_free_energy_per_link_0 = model.relative_gibbs_free_energy_per_link(&(ZERO*BOLTZMANN_CONSTANT*temperature/link_length), &temperature);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_relative_gibbs_free_energy_0 = model.nondimensional_relative_gibbs_free_energy(&ZERO);
            assert!(nondimensional_relative_gibbs_free_energy_0.abs() <= ZERO);
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
            let well_width = parameters.well_width_reference + parameters.well_width_scale*(0.5 - rng.gen::<f64>());
            let model = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
            let nondimensional_relative_gibbs_free_energy_per_link_0 = model.nondimensional_relative_gibbs_free_energy_per_link(&ZERO);
            assert!(nondimensional_relative_gibbs_free_energy_per_link_0.abs() <= ZERO);
        }
    }
}
