#![cfg(test)]
use super::*;
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
    pub nondimensional_potential_distance_reference: f64,
    pub nondimensional_potential_distance_scale: f64,
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
            abs_tol: 1e-8,
            rel_tol: 1e-6,
            number_of_loops: 88,
            hinge_mass_reference: 1e0,
            hinge_mass_scale: 1e0,
            link_length_reference: 1e0,
            link_length_scale: 1e0,
            number_of_links_minimum: 6,
            number_of_links_maximum: 25,
            nondimensional_potential_distance_reference: 1e0,
            nondimensional_potential_distance_scale: 2e0,
            nondimensional_potential_stiffness_reference: 5e2,
            nondimensional_potential_stiffness_scale: 1e3,
            temperature_reference: 3e2,
            temperature_scale: 1e2,
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
// mod nondimensional
// {
//     use super::*;
//     use rand::Rng;
//     #[test]
//     fn end_to_end_length()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
//             let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let end_to_end_length = model.end_to_end_length(&potential_distance, &potential_stiffness, &temperature);
//             let residual_abs = &end_to_end_length/link_length - &nondimensional_end_to_end_length;
//             let residual_rel = &residual_abs/&nondimensional_end_to_end_length;
//             assert!(residual_abs.abs() <= parameters.abs_tol);
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
//     #[test]
//     fn end_to_end_length_per_link()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
//             let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let end_to_end_length = model.end_to_end_length_per_link(&potential_distance, &potential_stiffness, &temperature);
//             let residual_abs = &end_to_end_length/link_length - &nondimensional_end_to_end_length;
//             let residual_rel = &residual_abs/&nondimensional_end_to_end_length;
//             assert!(residual_abs.abs() <= parameters.abs_tol);
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
//     #[test]
//     fn gibbs_free_energy()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
//             let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let gibbs_free_energy = model.gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature);
//             let residual_abs = gibbs_free_energy/BOLTZMANN_CONSTANT/temperature - nondimensional_gibbs_free_energy;
//             let residual_rel = residual_abs/nondimensional_gibbs_free_energy;
//             assert!(residual_abs.abs() <= parameters.abs_tol);
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
//     #[test]
//     fn gibbs_free_energy_per_link()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
//             let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
//             let residual_abs = gibbs_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - nondimensional_gibbs_free_energy_per_link;
//             let residual_rel = residual_abs/nondimensional_gibbs_free_energy_per_link;
//             assert!(residual_abs.abs() <= parameters.abs_tol);
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
//     #[test]
//     fn relative_gibbs_free_energy()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
//             let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature);
//             let residual_abs = relative_gibbs_free_energy/BOLTZMANN_CONSTANT/temperature - nondimensional_relative_gibbs_free_energy;
//             let residual_rel = residual_abs/nondimensional_relative_gibbs_free_energy;
//             assert!(residual_abs.abs() <= parameters.abs_tol);
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
//     #[test]
//     fn relative_gibbs_free_energy_per_link()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
//             let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
//             let residual_abs = relative_gibbs_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - nondimensional_relative_gibbs_free_energy_per_link;
//             let residual_rel = residual_abs/nondimensional_relative_gibbs_free_energy_per_link;
//             assert!(residual_abs.abs() <= parameters.abs_tol);
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
// }
// mod per_link
// {
//     use super::*;
//     use rand::Rng;
//     #[test]
//     fn end_to_end_length()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let end_to_end_length = model.end_to_end_length(&potential_distance, &potential_stiffness, &temperature);
//             let end_to_end_length_per_link = model.end_to_end_length_per_link(&potential_distance, &potential_stiffness, &temperature);
//             let residual_abs = end_to_end_length/(number_of_links as f64) - end_to_end_length_per_link;
//             let residual_rel = residual_abs/end_to_end_length_per_link;
//             assert!(residual_abs.abs() <= parameters.abs_tol);
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
//     #[test]
//     fn nondimensional_end_to_end_length()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
//             let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
//             let residual_abs = nondimensional_end_to_end_length/(number_of_links as f64) - nondimensional_end_to_end_length_per_link;
//             let residual_rel = residual_abs/nondimensional_end_to_end_length_per_link;
//             assert!(residual_abs.abs() <= parameters.abs_tol);
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
//     #[test]
//     fn gibbs_free_energy()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let gibbs_free_energy = model.gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature);
//             let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
//             let residual_abs = gibbs_free_energy/(number_of_links as f64) - gibbs_free_energy_per_link;
//             let residual_rel = residual_abs/gibbs_free_energy_per_link;
//             assert!(residual_abs.abs() <= parameters.abs_tol);
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
//     #[test]
//     fn relative_gibbs_free_energy()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature);
//             let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
//             let residual_abs = relative_gibbs_free_energy/(number_of_links as f64) - relative_gibbs_free_energy_per_link;
//             let residual_rel = residual_abs/relative_gibbs_free_energy_per_link;
//             assert!(residual_abs.abs() <= parameters.abs_tol);
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
//     #[test]
//     fn nondimensional_gibbs_free_energy()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
//             let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
//             let residual_abs = nondimensional_gibbs_free_energy/(number_of_links as f64) - nondimensional_gibbs_free_energy_per_link;
//             let residual_rel = residual_abs/nondimensional_gibbs_free_energy_per_link;
//             assert!(residual_abs.abs() <= parameters.abs_tol);
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
//     #[test]
//     fn nondimensional_relative_gibbs_free_energy()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
//             let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
//             let residual_abs = nondimensional_relative_gibbs_free_energy/(number_of_links as f64) - nondimensional_relative_gibbs_free_energy_per_link;
//             let residual_rel = residual_abs/nondimensional_relative_gibbs_free_energy_per_link;
//             assert!(residual_abs.abs() <= parameters.abs_tol);
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
// }
// mod relative
// {
//     use super::*;
//     use rand::Rng;
//     #[test]
//     fn gibbs_free_energy()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let gibbs_free_energy = model.gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature);
//             let gibbs_free_energy_0 = model.gibbs_free_energy(&ZERO, &potential_stiffness, &temperature);
//             let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&potential_distance, &potential_stiffness, &temperature);
//             let residual_abs = &gibbs_free_energy - &gibbs_free_energy_0 - &relative_gibbs_free_energy;
//             let residual_rel = &residual_abs/&gibbs_free_energy_0;
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
//     #[test]
//     fn gibbs_free_energy_per_link()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let potential_distance = nondimensional_potential_distance*(number_of_links as f64)*link_length;
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
//             let gibbs_free_energy_per_link_0 = model.gibbs_free_energy_per_link(&ZERO, &potential_stiffness, &temperature);
//             let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&potential_distance, &potential_stiffness, &temperature);
//             let residual_abs = &gibbs_free_energy_per_link - &gibbs_free_energy_per_link_0 - &relative_gibbs_free_energy_per_link;
//             let residual_rel = &residual_abs/&gibbs_free_energy_per_link_0;
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
//     #[test]
//     fn nondimensional_gibbs_free_energy()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
//             let nondimensional_gibbs_free_energy_0 = model.nondimensional_gibbs_free_energy(&ZERO, &nondimensional_potential_stiffness, &temperature);
//             let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
//             let residual_abs = &nondimensional_gibbs_free_energy - &nondimensional_gibbs_free_energy_0 - &nondimensional_relative_gibbs_free_energy;
//             let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_0;
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
//     #[test]
//     fn nondimensional_gibbs_free_energy_per_link()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_distance = parameters.nondimensional_potential_distance_reference + parameters.nondimensional_potential_distance_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness, &temperature);
//             let nondimensional_gibbs_free_energy_per_link_0 = model.nondimensional_gibbs_free_energy_per_link(&ZERO, &nondimensional_potential_stiffness, &temperature);
//             let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_potential_distance, &nondimensional_potential_stiffness);
//             let residual_abs = &nondimensional_gibbs_free_energy_per_link - &nondimensional_gibbs_free_energy_per_link_0 - &nondimensional_relative_gibbs_free_energy_per_link;
//             let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link_0;
//             assert!(residual_rel.abs() <= parameters.rel_tol);
//         }
//     }
// }
// mod zero
// {
//     use super::*;
//     use rand::Rng;
//     #[test]
//     fn end_to_end_length()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let end_to_end_length_0 = model.end_to_end_length(&ZERO, &potential_stiffness, &temperature);
//             assert!(end_to_end_length_0.abs() <= ZERO);
//         }
//     }
//     #[test]
//     fn end_to_end_length_per_link()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let end_to_end_length_0 = model.end_to_end_length_per_link(&ZERO, &potential_stiffness, &temperature);
//             assert!(end_to_end_length_0.abs() <= ZERO);
//         }
//     }
//     #[test]
//     fn nondimensional_end_to_end_length()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let end_to_end_length_0 = model.nondimensional_end_to_end_length(&ZERO, &nondimensional_potential_stiffness);
//             assert!(end_to_end_length_0.abs() <= ZERO);
//         }
//     }
//     #[test]
//     fn nondimensional_end_to_end_length_per_link()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let end_to_end_length_0 = model.nondimensional_end_to_end_length_per_link(&ZERO, &nondimensional_potential_stiffness);
//             assert!(end_to_end_length_0.abs() <= ZERO);
//         }
//     }
//     #[test]
//     fn relative_gibbs_free_energy()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let relative_gibbs_free_energy_0 = model.relative_gibbs_free_energy(&ZERO, &potential_stiffness, &temperature);
//             assert!(relative_gibbs_free_energy_0.abs() <= ZERO);
//         }
//     }
//     #[test]
//     fn relative_gibbs_free_energy_per_link()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
//             let potential_stiffness = nondimensional_potential_stiffness/((number_of_links as f64)*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature;
//             let relative_gibbs_free_energy_per_link_0 = model.relative_gibbs_free_energy_per_link(&ZERO, &potential_stiffness, &temperature);
//             assert!(relative_gibbs_free_energy_per_link_0.abs() <= ZERO);
//         }
//     }
//     #[test]
//     fn nondimensional_relative_gibbs_free_energy()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_relative_gibbs_free_energy_0 = model.nondimensional_relative_gibbs_free_energy(&ZERO, &nondimensional_potential_stiffness);
//             assert!(nondimensional_relative_gibbs_free_energy_0.abs() <= ZERO);
//         }
//     }
//     #[test]
//     fn nondimensional_relative_gibbs_free_energy_per_link()
//     {
//         let mut rng = rand::thread_rng();
//         let parameters = Parameters::default();
//         for _ in 0..parameters.number_of_loops
//         {
//             let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
//             let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
//             let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
//             let model = FJC::init(number_of_links, link_length, hinge_mass);
//             let nondimensional_potential_stiffness = parameters.nondimensional_potential_stiffness_reference + parameters.nondimensional_potential_stiffness_scale*(0.5 - rng.gen::<f64>());
//             let nondimensional_relative_gibbs_free_energy_per_link_0 = model.nondimensional_relative_gibbs_free_energy_per_link(&ZERO, &nondimensional_potential_stiffness);
//             assert!(nondimensional_relative_gibbs_free_energy_per_link_0.abs() <= ZERO);
//         }
//     }
// }