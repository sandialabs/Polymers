#![cfg(test)]
use super::{
    ideal::Ideal,
    fjc::FJC,
    efjc::EFJC,
    swfjc::SWFJC
};
use crate::physics::BOLTZMANN_CONSTANT;
pub fn integrate<F>(function: F, lower_lim: &f64, upper_lim: &f64, num_points: &u128) -> f64
where F: Fn(f64) -> f64
{
    let dx = (upper_lim - lower_lim)/(*num_points as f64);
    (0..=num_points-1).collect::<Vec::<u128>>().iter().map(|index| function(lower_lim + (0.5 + *index as f64)*dx)).sum::<f64>()*dx
}
pub struct Parameters
{
    pub abs_tol: f64,
    pub rel_tol: f64,
    pub rel_tol_thermodynamic_limit: f64,
    pub log_log_tol: f64,
    pub log_log_scale: f64,
    pub number_of_loops: u32,
    pub hinge_mass_reference: f64,
    pub hinge_mass_scale: f64,
    pub link_length_reference: f64,
    pub link_length_scale: f64,
    pub number_of_links_minimum: u8,
    pub number_of_links_maximum: u8,
    pub link_stiffness_reference: f64,
    pub link_stiffness_scale: f64,
    pub nondimensional_link_stiffness_large: f64,
    pub nondimensional_link_stiffness_medium: f64,
    pub well_width_reference: f64,
    pub well_width_scale: f64,
    pub nondimensional_well_width_small: f64,
    pub nondimensional_end_to_end_length_per_link_reference: f64,
    pub nondimensional_end_to_end_length_per_link_scale: f64,
    pub nondimensional_end_to_end_length_per_link_small: f64,
    pub nondimensional_force_reference: f64,
    pub nondimensional_force_scale: f64,
    pub nondimensional_force_small: f64,
    pub nondimensional_potential_distance_reference: f64,
    pub nondimensional_potential_distance_scale: f64,
    pub nondimensional_potential_distance_small: f64,
    pub nondimensional_potential_distance_large_1: f64,
    pub nondimensional_potential_distance_large_2: f64,
    pub nondimensional_potential_stiffness_reference: f64,
    pub nondimensional_potential_stiffness_scale: f64,
    pub nondimensional_potential_stiffness_small: f64,
    pub nondimensional_potential_stiffness_large: f64,
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
            rel_tol_thermodynamic_limit: 1e-1,
            log_log_tol: 5e-2,
            log_log_scale: 12e-1,
            number_of_loops: 8,
            hinge_mass_reference: 1e0,
            hinge_mass_scale: 1e0,
            link_length_reference: 1e0,
            link_length_scale: 1e0,
            number_of_links_minimum: 5,
            number_of_links_maximum: 25,
            link_stiffness_reference: 5e5,
            link_stiffness_scale: 99e4,
            nondimensional_link_stiffness_large: 1e4,
            nondimensional_link_stiffness_medium: 1e1,
            well_width_reference: 99e-2,
            well_width_scale: 5e-1,
            nondimensional_well_width_small: 1e-2,
            nondimensional_end_to_end_length_per_link_reference: 5e-1,
            nondimensional_end_to_end_length_per_link_scale: 99e-2,
            nondimensional_end_to_end_length_per_link_small: 25e-2,
            nondimensional_force_reference: 5e1,
            nondimensional_force_scale: 1e2,
            nondimensional_force_small: 75e-2,
            nondimensional_potential_distance_reference: 1e0,
            nondimensional_potential_distance_scale: 2e0,
            nondimensional_potential_distance_small: 25e-2,
            nondimensional_potential_distance_large_1: 1e1,
            nondimensional_potential_distance_large_2: 1e1 + 25e-1,
            nondimensional_potential_stiffness_reference: 5e1,
            nondimensional_potential_stiffness_scale: 1e2,
            nondimensional_potential_stiffness_small: 1e-2,
            nondimensional_potential_stiffness_large: 1e2,
            temperature_reference: 3e2,
            temperature_scale: 1e2,
        }
    }
}
mod fjc
{
    use super::*;
    mod ideal
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let force = BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force;
                let end_to_end_length_fjc = fjc.thermodynamics.isotensional.end_to_end_length(&force, &temperature);
                let end_to_end_length_ideal = ideal.thermodynamics.isotensional.end_to_end_length(&force, &temperature);
                let residual_abs = &end_to_end_length_fjc - &end_to_end_length_ideal;
                let residual_rel = &residual_abs/&end_to_end_length_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let force = BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force;
                let end_to_end_length_per_link_fjc = fjc.thermodynamics.isotensional.end_to_end_length_per_link(&force, &temperature);
                let end_to_end_length_per_link_ideal = ideal.thermodynamics.isotensional.end_to_end_length_per_link(&force, &temperature);
                let residual_abs = &end_to_end_length_per_link_fjc - &end_to_end_length_per_link_ideal;
                let residual_rel = &residual_abs/&end_to_end_length_per_link_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let nondimensional_end_to_end_length_fjc = fjc.thermodynamics.isotensional.nondimensional_end_to_end_length(&nondimensional_force);
                let nondimensional_end_to_end_length_ideal = ideal.thermodynamics.isotensional.nondimensional_end_to_end_length(&nondimensional_force);
                let residual_abs = &nondimensional_end_to_end_length_fjc - &nondimensional_end_to_end_length_ideal;
                let residual_rel = &residual_abs/&nondimensional_end_to_end_length_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let nondimensional_end_to_end_length_per_link_fjc = fjc.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                let nondimensional_end_to_end_length_per_link_ideal = ideal.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                let residual_abs = &nondimensional_end_to_end_length_per_link_fjc - &nondimensional_end_to_end_length_per_link_ideal;
                let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
                let force_fjc = fjc.thermodynamics.isometric.force(&end_to_end_length, &temperature);
                let force_ideal = ideal.thermodynamics.isometric.force(&end_to_end_length, &temperature);
                let residual_abs = &force_fjc - &force_ideal;
                let residual_rel = &residual_abs/&force_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_end_to_end_length_per_link);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_small*(1.0 - 0.5*rng.gen::<f64>());
                let nondimensional_force_fjc = fjc.thermodynamics.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link);
                let nondimensional_force_ideal = ideal.thermodynamics.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link);
                let residual_abs = &nondimensional_force_fjc - &nondimensional_force_ideal;
                let residual_rel = &residual_abs/&nondimensional_force_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_end_to_end_length_per_link);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
                let helmholtz_free_energy_fjc = fjc.thermodynamics.isometric.helmholtz_free_energy(&end_to_end_length, &temperature);
                let helmholtz_free_energy_ideal = ideal.thermodynamics.isometric.helmholtz_free_energy(&end_to_end_length, &temperature);
                let residual_abs = &helmholtz_free_energy_fjc - &helmholtz_free_energy_ideal;
                let residual_rel = &residual_abs/&helmholtz_free_energy_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_end_to_end_length_per_link);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
                let helmholtz_free_energy_per_link_fjc = fjc.thermodynamics.isometric.helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
                let helmholtz_free_energy_per_link_ideal = ideal.thermodynamics.isometric.helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
                let residual_abs = &helmholtz_free_energy_per_link_fjc - &helmholtz_free_energy_per_link_ideal;
                let residual_rel = &residual_abs/&helmholtz_free_energy_per_link_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_end_to_end_length_per_link);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
                let relative_helmholtz_free_energy_fjc = fjc.thermodynamics.isometric.relative_helmholtz_free_energy(&end_to_end_length, &temperature);
                let relative_helmholtz_free_energy_ideal = ideal.thermodynamics.isometric.relative_helmholtz_free_energy(&end_to_end_length, &temperature);
                let residual_abs = &relative_helmholtz_free_energy_fjc - &relative_helmholtz_free_energy_ideal;
                let residual_rel = &residual_abs/&relative_helmholtz_free_energy_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_end_to_end_length_per_link);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
                let relative_helmholtz_free_energy_per_link_fjc = fjc.thermodynamics.isometric.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
                let relative_helmholtz_free_energy_per_link_ideal = ideal.thermodynamics.isometric.relative_helmholtz_free_energy_per_link(&end_to_end_length, &temperature);
                let residual_abs = &relative_helmholtz_free_energy_per_link_fjc - &relative_helmholtz_free_energy_per_link_ideal;
                let residual_rel = &residual_abs/&relative_helmholtz_free_energy_per_link_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_end_to_end_length_per_link);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let nondimensional_helmholtz_free_energy_fjc = fjc.thermodynamics.isometric.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
                let nondimensional_helmholtz_free_energy_ideal = ideal.thermodynamics.isometric.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, &temperature);
                let residual_abs = &nondimensional_helmholtz_free_energy_fjc - &nondimensional_helmholtz_free_energy_ideal;
                let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_end_to_end_length_per_link);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let nondimensional_helmholtz_free_energy_per_link_fjc = fjc.thermodynamics.isometric.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
                let nondimensional_helmholtz_free_energy_per_link_ideal = ideal.thermodynamics.isometric.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, &temperature);
                let residual_abs = &nondimensional_helmholtz_free_energy_per_link_fjc - &nondimensional_helmholtz_free_energy_per_link_ideal;
                let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_end_to_end_length_per_link);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_small*(1.0 - 0.5*rng.gen::<f64>());
                let nondimensional_relative_helmholtz_free_energy_fjc = fjc.thermodynamics.isometric.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link);
                let nondimensional_relative_helmholtz_free_energy_ideal = ideal.thermodynamics.isometric.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link);
                let residual_abs = &nondimensional_relative_helmholtz_free_energy_fjc - &nondimensional_relative_helmholtz_free_energy_ideal;
                let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_end_to_end_length_per_link);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_small*(1.0 - 0.5*rng.gen::<f64>());
                let nondimensional_relative_helmholtz_free_energy_per_link_fjc = fjc.thermodynamics.isometric.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
                let nondimensional_relative_helmholtz_free_energy_per_link_ideal = ideal.thermodynamics.isometric.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
                let residual_abs = &nondimensional_relative_helmholtz_free_energy_per_link_fjc - &nondimensional_relative_helmholtz_free_energy_per_link_ideal;
                let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_end_to_end_length_per_link);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let force = BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force;
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let gibbs_free_energy_fjc = fjc.thermodynamics.isotensional.gibbs_free_energy(&force, &temperature);
                let gibbs_free_energy_ideal = ideal.thermodynamics.isotensional.gibbs_free_energy(&force, &temperature);
                let residual_abs = &gibbs_free_energy_fjc - &gibbs_free_energy_ideal;
                let residual_rel = &residual_abs/&gibbs_free_energy_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let force = BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force;
                let gibbs_free_energy_per_link_fjc = fjc.thermodynamics.isotensional.gibbs_free_energy_per_link(&force, &temperature);
                let gibbs_free_energy_per_link_ideal = ideal.thermodynamics.isotensional.gibbs_free_energy_per_link(&force, &temperature);
                let residual_abs = &gibbs_free_energy_per_link_fjc - &gibbs_free_energy_per_link_ideal;
                let residual_rel = &residual_abs/&gibbs_free_energy_per_link_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let force = BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force;
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let relative_gibbs_free_energy_fjc = fjc.thermodynamics.isotensional.relative_gibbs_free_energy(&force, &temperature);
                let relative_gibbs_free_energy_ideal = ideal.thermodynamics.isotensional.relative_gibbs_free_energy(&force, &temperature);
                let residual_abs = &relative_gibbs_free_energy_fjc - &relative_gibbs_free_energy_ideal;
                let residual_rel = &residual_abs/&relative_gibbs_free_energy_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let force = BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force;
                let relative_gibbs_free_energy_per_link_fjc = fjc.thermodynamics.isotensional.relative_gibbs_free_energy_per_link(&force, &temperature);
                let relative_gibbs_free_energy_per_link_ideal = ideal.thermodynamics.isotensional.relative_gibbs_free_energy_per_link(&force, &temperature);
                let residual_abs = &relative_gibbs_free_energy_per_link_fjc - &relative_gibbs_free_energy_per_link_ideal;
                let residual_rel = &residual_abs/&relative_gibbs_free_energy_per_link_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let nondimensional_gibbs_free_energy_fjc = fjc.thermodynamics.isotensional.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature);
                let nondimensional_gibbs_free_energy_ideal = ideal.thermodynamics.isotensional.nondimensional_gibbs_free_energy(&nondimensional_force, &temperature);
                let residual_abs = &nondimensional_gibbs_free_energy_fjc - &nondimensional_gibbs_free_energy_ideal;
                let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let nondimensional_gibbs_free_energy_per_link_fjc = fjc.thermodynamics.isotensional.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature);
                let nondimensional_gibbs_free_energy_per_link_ideal = ideal.thermodynamics.isotensional.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, &temperature);
                let residual_abs = &nondimensional_gibbs_free_energy_per_link_fjc - &nondimensional_gibbs_free_energy_per_link_ideal;
                let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let nondimensional_relative_gibbs_free_energy_fjc = fjc.thermodynamics.isotensional.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
                let nondimensional_relative_gibbs_free_energy_ideal = ideal.thermodynamics.isotensional.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
                let residual_abs = &nondimensional_relative_gibbs_free_energy_fjc - &nondimensional_relative_gibbs_free_energy_ideal;
                let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let nondimensional_relative_gibbs_free_energy_per_link_fjc = fjc.thermodynamics.isotensional.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
                let nondimensional_relative_gibbs_free_energy_per_link_ideal = ideal.thermodynamics.isotensional.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
                let residual_abs = &nondimensional_relative_gibbs_free_energy_per_link_fjc - &nondimensional_relative_gibbs_free_energy_per_link_ideal;
                let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
            }
        }
        #[test]
        fn equilibrium_distribution()
        {
            let mut rng = rand::thread_rng();
            let parameters = Parameters::default();
            for _ in 0..parameters.number_of_loops
            {
                let number_of_links: u8 = parameters.number_of_links_maximum;
                let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_small*(1.0 - 0.5*rng.gen::<f64>());
                let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
                let equilibrium_distribution_fjc = fjc.thermodynamics.isometric.equilibrium_distribution(&end_to_end_length);
                let equilibrium_distribution_ideal = ideal.thermodynamics.isometric.equilibrium_distribution(&end_to_end_length);
                let residual_abs = &equilibrium_distribution_fjc - &equilibrium_distribution_ideal;
                let residual_rel = &residual_abs/&equilibrium_distribution_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_end_to_end_length_per_link);
            }
        }
        #[test]
        fn nondimensional_equilibrium_distribution()
        {
            let mut rng = rand::thread_rng();
            let parameters = Parameters::default();
            for _ in 0..parameters.number_of_loops
            {
                let number_of_links: u8 = parameters.number_of_links_maximum;
                let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_small*(1.0 - 0.5*rng.gen::<f64>());
                let nondimensional_equilibrium_distribution_fjc = fjc.thermodynamics.isometric.nondimensional_equilibrium_distribution(&nondimensional_end_to_end_length_per_link);
                let nondimensional_equilibrium_distribution_ideal = ideal.thermodynamics.isometric.nondimensional_equilibrium_distribution(&nondimensional_end_to_end_length_per_link);
                let residual_abs = &nondimensional_equilibrium_distribution_fjc - &nondimensional_equilibrium_distribution_ideal;
                let residual_rel = &residual_abs/&nondimensional_equilibrium_distribution_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_end_to_end_length_per_link);
            }
        }
        #[test]
        fn equilibrium_radial_distribution()
        {
            let mut rng = rand::thread_rng();
            let parameters = Parameters::default();
            for _ in 0..parameters.number_of_loops
            {
                let number_of_links: u8 = parameters.number_of_links_maximum;
                let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_small*(1.0 - 0.5*rng.gen::<f64>());
                let end_to_end_length = nondimensional_end_to_end_length_per_link*(number_of_links as f64)*link_length;
                let equilibrium_radial_distribution_fjc = fjc.thermodynamics.isometric.equilibrium_radial_distribution(&end_to_end_length);
                let equilibrium_radial_distribution_ideal = ideal.thermodynamics.isometric.equilibrium_radial_distribution(&end_to_end_length);
                let residual_abs = &equilibrium_radial_distribution_fjc - &equilibrium_radial_distribution_ideal;
                let residual_rel = &residual_abs/&equilibrium_radial_distribution_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_end_to_end_length_per_link);
            }
        }
        #[test]
        fn nondimensional_equilibrium_radial_distribution()
        {
            let mut rng = rand::thread_rng();
            let parameters = Parameters::default();
            for _ in 0..parameters.number_of_loops
            {
                let number_of_links: u8 = parameters.number_of_links_maximum;
                let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_small*(1.0 - 0.5*rng.gen::<f64>());
                let nondimensional_equilibrium_radial_distribution_fjc = fjc.thermodynamics.isometric.nondimensional_equilibrium_radial_distribution(&nondimensional_end_to_end_length_per_link);
                let nondimensional_equilibrium_radial_distribution_ideal = ideal.thermodynamics.isometric.nondimensional_equilibrium_radial_distribution(&nondimensional_end_to_end_length_per_link);
                let residual_abs = &nondimensional_equilibrium_radial_distribution_fjc - &nondimensional_equilibrium_radial_distribution_ideal;
                let residual_rel = &residual_abs/&nondimensional_equilibrium_radial_distribution_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_end_to_end_length_per_link);
            }
        }
    }
}
mod efjc
{
    use super::*;
    mod ideal
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
                let link_stiffness = parameters.link_stiffness_scale;
                let efjc = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let force = BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force;
                let end_to_end_length_efjc = efjc.thermodynamics.isotensional.end_to_end_length(&force, &temperature);
                let end_to_end_length_ideal = ideal.thermodynamics.isotensional.end_to_end_length(&force, &temperature);
                let residual_abs = &end_to_end_length_efjc - &end_to_end_length_ideal;
                let residual_rel = &residual_abs/&end_to_end_length_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let link_stiffness = parameters.link_stiffness_scale;
                let efjc = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let force = BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force;
                let end_to_end_length_per_link_efjc = efjc.thermodynamics.isotensional.end_to_end_length_per_link(&force, &temperature);
                let end_to_end_length_per_link_ideal = ideal.thermodynamics.isotensional.end_to_end_length_per_link(&force, &temperature);
                let residual_abs = &end_to_end_length_per_link_efjc - &end_to_end_length_per_link_ideal;
                let residual_rel = &residual_abs/&end_to_end_length_per_link_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let link_stiffness = parameters.link_stiffness_scale;
                let efjc = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let nondimensional_end_to_end_length_efjc = efjc.thermodynamics.isotensional.nondimensional_end_to_end_length(&nondimensional_force, &temperature);
                let nondimensional_end_to_end_length_ideal = ideal.thermodynamics.isotensional.nondimensional_end_to_end_length(&nondimensional_force);
                let residual_abs = &nondimensional_end_to_end_length_efjc - &nondimensional_end_to_end_length_ideal;
                let residual_rel = &residual_abs/&nondimensional_end_to_end_length_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let link_stiffness = parameters.link_stiffness_scale;
                let efjc = EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let nondimensional_end_to_end_length_per_link_efjc = efjc.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force, &temperature);
                let nondimensional_end_to_end_length_per_link_ideal = ideal.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                let residual_abs = &nondimensional_end_to_end_length_per_link_efjc - &nondimensional_end_to_end_length_per_link_ideal;
                let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
            }
        }
    }
}
mod swfjc
{
    use super::*;
    mod ideal
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
                let well_width = parameters.nondimensional_well_width_small*link_length;
                let swfjc = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let force = BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force;
                let end_to_end_length_swfjc = swfjc.thermodynamics.isotensional.end_to_end_length(&force, &temperature);
                let end_to_end_length_ideal = ideal.thermodynamics.isotensional.end_to_end_length(&force, &temperature);
                let residual_abs = &end_to_end_length_swfjc - &end_to_end_length_ideal;
                let residual_rel = &residual_abs/&end_to_end_length_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let well_width = parameters.nondimensional_well_width_small*link_length;
                let swfjc = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let force = BOLTZMANN_CONSTANT*temperature/link_length*nondimensional_force;
                let end_to_end_length_per_link_swfjc = swfjc.thermodynamics.isotensional.end_to_end_length_per_link(&force, &temperature);
                let end_to_end_length_per_link_ideal = ideal.thermodynamics.isotensional.end_to_end_length_per_link(&force, &temperature);
                let residual_abs = &end_to_end_length_per_link_swfjc - &end_to_end_length_per_link_ideal;
                let residual_rel = &residual_abs/&end_to_end_length_per_link_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let well_width = parameters.nondimensional_well_width_small*link_length;
                let swfjc = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let nondimensional_end_to_end_length_swfjc = swfjc.thermodynamics.isotensional.nondimensional_end_to_end_length(&nondimensional_force);
                let nondimensional_end_to_end_length_ideal = ideal.thermodynamics.isotensional.nondimensional_end_to_end_length(&nondimensional_force);
                let residual_abs = &nondimensional_end_to_end_length_swfjc - &nondimensional_end_to_end_length_ideal;
                let residual_rel = &residual_abs/&nondimensional_end_to_end_length_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
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
                let well_width = parameters.nondimensional_well_width_small*link_length;
                let swfjc = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
                let ideal = Ideal::init(number_of_links, link_length, hinge_mass);
                let nondimensional_force = parameters.nondimensional_force_small*(1.0 - 0.5*rng.gen::<f64>());
                let nondimensional_end_to_end_length_per_link_swfjc = swfjc.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                let nondimensional_end_to_end_length_per_link_ideal = ideal.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                let residual_abs = &nondimensional_end_to_end_length_per_link_swfjc - &nondimensional_end_to_end_length_per_link_ideal;
                let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link_ideal;
                assert!(residual_rel.abs() <= parameters.rel_tol_thermodynamic_limit);
                assert!(residual_rel.abs() <= nondimensional_force);
            }
        }
    }
    mod fjc
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
                let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                let link_length = rng.gen::<f64>();
                let well_width = parameters.nondimensional_well_width_small*link_length;
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let swfjc = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
                let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
                let end_to_end_length_fjc = fjc.thermodynamics.isotensional.end_to_end_length(&force, &temperature);
                let end_to_end_length_swfjc = swfjc.thermodynamics.isotensional.end_to_end_length(&force, &temperature);
                let residual = &end_to_end_length_swfjc - &end_to_end_length_fjc;
                assert!(residual.abs() <= (number_of_links as f64)*link_length*parameters.nondimensional_well_width_small);
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
                let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                let link_length = rng.gen::<f64>();
                let well_width = parameters.nondimensional_well_width_small*link_length;
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let swfjc = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
                let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                let force = nondimensional_force*BOLTZMANN_CONSTANT*temperature/link_length;
                let end_to_end_length_per_link_fjc = fjc.thermodynamics.isotensional.end_to_end_length_per_link(&force, &temperature);
                let end_to_end_length_per_link_swfjc = swfjc.thermodynamics.isotensional.end_to_end_length_per_link(&force, &temperature);
                let residual = &end_to_end_length_per_link_swfjc - &end_to_end_length_per_link_fjc;
                assert!(residual.abs() <= link_length*parameters.nondimensional_well_width_small);
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
                let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                let link_length = rng.gen::<f64>();
                let well_width = parameters.nondimensional_well_width_small*link_length;
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let swfjc = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
                let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                let nondimensional_end_to_end_length_fjc = fjc.thermodynamics.isotensional.nondimensional_end_to_end_length(&nondimensional_force);
                let nondimensional_end_to_end_length_swfjc = swfjc.thermodynamics.isotensional.nondimensional_end_to_end_length(&nondimensional_force);
                let residual = &nondimensional_end_to_end_length_swfjc - &nondimensional_end_to_end_length_fjc;
                assert!(residual.abs() <= (number_of_links as f64)*parameters.nondimensional_well_width_small);
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
                let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                let link_length = rng.gen::<f64>();
                let well_width = parameters.nondimensional_well_width_small*link_length;
                let fjc = FJC::init(number_of_links, link_length, hinge_mass);
                let swfjc = SWFJC::init(number_of_links, link_length, hinge_mass, well_width);
                let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                let nondimensional_end_to_end_length_per_link_fjc = fjc.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                let nondimensional_end_to_end_length_per_link_swfjc = swfjc.thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                let residual = &nondimensional_end_to_end_length_per_link_swfjc - &nondimensional_end_to_end_length_per_link_fjc;
                assert!(residual.abs() <= parameters.nondimensional_well_width_small);
            }
        }
    }
}
