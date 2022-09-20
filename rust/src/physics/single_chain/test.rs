#![cfg(test)]

use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};

pub struct Parameters
{
    pub number_of_loops: u32,
    pub abs_tol_for_close: f64,
    pub rel_tol_for_close: f64,
    pub abs_tol_for_equals: f64,
    pub rel_tol_for_equals: f64,
    pub force_reference: f64,
    pub force_scale: f64,
    pub nondimensional_force_reference: f64,
    pub nondimensional_force_scale: f64,
    pub temperature_reference: f64,
    pub temperature_scale: f64,
    pub hinge_mass_reference: f64,
    pub hinge_mass_scale: f64,
    pub link_length_reference: f64,
    pub link_length_scale: f64,
    pub default_link_length: f64,
    pub minimum_number_of_links: u16,
    pub maximum_number_of_links: u16,
    pub default_number_of_links: u16,
    pub default_hinge_mass: f64,
    pub zero: f64,
    pub small: f64
}

impl Default for Parameters
{
    fn default() -> Self
    {
        Self
        {
            number_of_loops: 888888,
            abs_tol_for_close: 1e-5,
            rel_tol_for_close: 1e-4,
            abs_tol_for_equals: 1e-9,
            rel_tol_for_equals: 1e-9,
            force_reference: 1e0,
            force_scale: 1e0,
            nondimensional_force_reference: 1e0,
            nondimensional_force_scale: 1e0,
            temperature_reference: 3e2,
            temperature_scale: 1e0,
            hinge_mass_reference: 1e0,
            hinge_mass_scale: 1e0,
            link_length_reference: 1e0,
            link_length_scale: 1e0,
            default_link_length: 1e0,
            minimum_number_of_links: 8,
            maximum_number_of_links: 88,
            default_number_of_links: 8,
            default_hinge_mass: 1e0,
            zero: 1e-88,
            small: 1e-2
        }
    }
}

macro_rules! base
{
    ( $model:ty ) =>
    {
        use rand::Rng;
        use crate::physics::single_chain::test::Parameters;
        #[test]
        fn init()
        {
            let parameters = Parameters::default();
            let _ = <$model>::init(parameters.default_number_of_links, parameters.default_link_length, parameters.default_hinge_mass);
        }
        #[test]
        fn number_of_links()
        {
            let parameters = Parameters::default();
            let mut rng = rand::thread_rng();
            for _ in 0..parameters.number_of_loops
            {
                let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                assert_eq!(number_of_links, <$model>::init(number_of_links, parameters.default_link_length, parameters.default_hinge_mass).number_of_links);
            }
        }
        #[test]
        fn link_length()
        {
            let parameters = Parameters::default();
            let mut rng = rand::thread_rng();
            for _ in 0..parameters.number_of_loops
            {
                let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                assert_eq!(link_length, <$model>::init(parameters.default_number_of_links, link_length, parameters.default_hinge_mass).link_length);
            }
        }
        #[test]
        fn hinge_mass()
        {
            let parameters = Parameters::default();
            let mut rng = rand::thread_rng();
            for _ in 0..parameters.number_of_loops
            {
                let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                assert_eq!(hinge_mass, <$model>::init(parameters.default_number_of_links, parameters.default_link_length, hinge_mass).hinge_mass);
            }
        }
        #[test]
        fn number_of_links_and_link_length_and_hinge_mass()
        {
            let parameters = Parameters::default();
            let mut rng = rand::thread_rng();
            for _ in 0..parameters.number_of_loops
            {
                let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                let link_length = rng.gen::<f64>();
                assert_eq!(link_length, <$model>::init(number_of_links, link_length, hinge_mass).link_length);
            }
        }
    }
}
pub(crate) use base;

macro_rules! single_chain
{
    ( $model:ty ) =>
    {
        mod init
        {
            use super::*;
            use crate::physics::single_chain::test::Parameters;
            #[test]
            fn mechanics()
            {
                let parameters = Parameters::default();
                let _ = <$model>::init(parameters.default_number_of_links, parameters.default_link_length, parameters.default_hinge_mass).mechanics;
            }
            #[test]
            fn thermodynamics()
            {
                let parameters = Parameters::default();
                let _ = <$model>::init(parameters.default_number_of_links, parameters.default_link_length, parameters.default_hinge_mass).thermodynamics;
            }
        }
    }	
}
pub(crate) use single_chain;

macro_rules! mechanics
{
    ( $model:ty ) =>
    {
        mod init
        {
            use super::*;
            use crate::physics::single_chain::test::Parameters;
            #[test]
            fn number_of_configurations()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    assert_eq!(number_of_links + 1, model.random_configuration().len() as u16);
                }
            }
        }
    }
}
pub(crate) use mechanics;

macro_rules! thermodynamics
{
    ( $model:ty ) =>
    {
        mod init
        {
            use super::*;
            use crate::physics::single_chain::test::Parameters;
            #[test]
            fn isometric()
            {
                let parameters = Parameters::default();
                let _ = <$model>::init(parameters.default_number_of_links, parameters.default_link_length, parameters.default_hinge_mass).isometric;
            }
            #[test]
            fn isotensional()
            {
                let parameters = Parameters::default();
                let _ = <$model>::init(parameters.default_number_of_links, parameters.default_link_length, parameters.default_hinge_mass).isotensional;
            }
        }
        // mod legendre
        // {
        //     #[test]
        //     fn nondimensional_end_to_end_length()
        //     {}
        //     #[test]
        //     fn nondimensional_end_to_end_length_per_link()
        //     {}
        //     #[test]
        //     fn nondimensional_force()
        //     {}
        //     #[test]
        //     fn nondimensional_force_per_link()
        //     {}
        //     #[test]
        //     fn nondimensional_relative_helmholtz_free_energy()
        //     {}
        //     #[test]
        //     fn nondimensional_relative_helmholtz_free_energy_per_link()
        //     {}
        //     #[test]
        //     fn nondimensional_relative_gibbs_free_energy()
        //     {}
        //     #[test]
        //     fn nondimensional_relative_gibbs_free_energy_per_link()
        //     {}
        // }
    }
}
pub(crate) use thermodynamics;

macro_rules! isometric
{
    ( $model:ty ) =>
    {
        mod init
        {
            use super::*;
            use crate::physics::single_chain::test::Parameters;
            #[test]
            fn legendre()
            {
                let parameters = Parameters::default();
                let _ = <$model>::init(parameters.default_number_of_links, parameters.default_link_length, parameters.default_hinge_mass).legendre;
            }
        }
    }
}
pub(crate) use isometric;

macro_rules! isotensional
{
    ( $model:ty ) =>
    {
        mod init
        {
            use super::*;
            use crate::physics::single_chain::test::Parameters;
            #[test]
            fn legendre()
            {
                let parameters = Parameters::default();
                let _ = <$model>::init(parameters.default_number_of_links, parameters.default_link_length, parameters.default_hinge_mass).legendre;
            }
        }
        mod nondimensional
        {
            use super::*;
            use rand::Rng;
            use crate::physics::single_chain::test::Parameters;
            #[test]
            fn end_to_end_length()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = model.end_to_end_length(&force, temperature);
                    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                    let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_force);
                    let residual_abs = &end_to_end_length/link_length - &nondimensional_end_to_end_length;
                    let residual_rel = &residual_abs/&nondimensional_end_to_end_length;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn end_to_end_length_per_link()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length_per_link = model.end_to_end_length_per_link(&force, temperature);
                    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                    let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                    let residual_abs = &end_to_end_length_per_link/link_length - &nondimensional_end_to_end_length_per_link;
                    let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn gibbs_free_energy()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let gibbs_free_energy = model.gibbs_free_energy(&force, temperature);
                    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                    let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_force, temperature);
                    let residual_abs = &gibbs_free_energy/BOLTZMANN_CONSTANT/temperature - &nondimensional_gibbs_free_energy;
                    let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn relative_gibbs_free_energy()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&force, temperature);
                    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                    let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
                    let residual_abs = &relative_gibbs_free_energy/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_gibbs_free_energy;
                    let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn gibbs_free_energy_per_link()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&force, temperature);
                    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                    let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, temperature);
                    let residual_abs = &gibbs_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn relative_gibbs_free_energy_per_link()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&force, temperature);
                    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                    let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
                    let residual_abs = &relative_gibbs_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            mod legendre
            {
                use super::*;
                use rand::Rng;
                use crate::physics::single_chain::test::Parameters;
                #[test]
                fn helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let helmholtz_free_energy = model.legendre.helmholtz_free_energy(&force, temperature);
                        let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                        let nondimensional_helmholtz_free_energy = model.legendre.nondimensional_helmholtz_free_energy(&nondimensional_force, temperature);
                        let residual_abs = &helmholtz_free_energy/BOLTZMANN_CONSTANT/temperature - &nondimensional_helmholtz_free_energy;
                        let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                    }
                }
                #[test]
                fn relative_helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let relative_helmholtz_free_energy = model.legendre.relative_helmholtz_free_energy(&force, temperature);
                        let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                        let nondimensional_relative_helmholtz_free_energy = model.legendre.nondimensional_relative_helmholtz_free_energy(&nondimensional_force);
                        let residual_abs = &relative_helmholtz_free_energy/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_helmholtz_free_energy;
                        let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                    }
                }
                #[test]
                fn helmholtz_free_energy_per_link()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let helmholtz_free_energy_per_link = model.legendre.helmholtz_free_energy_per_link(&force, temperature);
                        let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                        let nondimensional_helmholtz_free_energy_per_link = model.legendre.nondimensional_helmholtz_free_energy_per_link(&nondimensional_force, temperature);
                        let residual_abs = &helmholtz_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_helmholtz_free_energy_per_link;
                        let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                    }
                }
                #[test]
                fn relative_helmholtz_free_energy_per_link()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let relative_helmholtz_free_energy_per_link = model.legendre.relative_helmholtz_free_energy_per_link(&force, temperature);
                        let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                        let nondimensional_relative_helmholtz_free_energy_per_link = model.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force);
                        let residual_abs = &relative_helmholtz_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_helmholtz_free_energy_per_link;
                        let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                    }
                }
            }
        }
        mod per_link
        {
            use super::*;
            use rand::Rng;
            use crate::physics::single_chain::test::Parameters;
            #[test]
            fn end_to_end_length()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = model.end_to_end_length(&force, temperature);
                    let end_to_end_length_per_link = model.end_to_end_length_per_link(&force, temperature);
                    let residual_abs = &end_to_end_length/(model.number_of_links as f64) - &end_to_end_length_per_link;
                    let residual_rel = &residual_abs/&end_to_end_length_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn nondimensional_end_to_end_length()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_force);
                    let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                    let residual_abs = &nondimensional_end_to_end_length/(model.number_of_links as f64) - &nondimensional_end_to_end_length_per_link;
                    let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn gibbs_free_energy()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let gibbs_free_energy = model.gibbs_free_energy(&force, temperature);
                    let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&force, temperature);
                    let residual_abs = &gibbs_free_energy/(model.number_of_links as f64) - &gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn relative_gibbs_free_energy()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&force, temperature);
                    let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&force, temperature);
                    let residual_abs = &relative_gibbs_free_energy/(model.number_of_links as f64) - &relative_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&relative_gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn nondimensional_gibbs_free_energy()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_force, temperature);
                    let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, temperature);
                    let residual_abs = &nondimensional_gibbs_free_energy/(model.number_of_links as f64) - &nondimensional_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn nondimensional_relative_gibbs_free_energy()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
                    let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
                    let residual_abs = &nondimensional_relative_gibbs_free_energy/(model.number_of_links as f64) - &nondimensional_relative_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            mod legendre
            {
                use super::*;
                use rand::Rng;
                use crate::physics::single_chain::test::Parameters;
                #[test]
                fn helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let helmholtz_free_energy = model.legendre.helmholtz_free_energy(&force, temperature);
                        let helmholtz_free_energy_per_link = model.legendre.helmholtz_free_energy_per_link(&force, temperature);
                        let residual_abs = &helmholtz_free_energy/(model.number_of_links as f64) - &helmholtz_free_energy_per_link;
                        let residual_rel = &residual_abs/&helmholtz_free_energy_per_link;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                    }
                }
                #[test]
                fn relative_helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let relative_helmholtz_free_energy = model.legendre.relative_helmholtz_free_energy(&force, temperature);
                        let relative_helmholtz_free_energy_per_link = model.legendre.relative_helmholtz_free_energy_per_link(&force, temperature);
                        let residual_abs = &relative_helmholtz_free_energy/(model.number_of_links as f64) - &relative_helmholtz_free_energy_per_link;
                        let residual_rel = &residual_abs/&relative_helmholtz_free_energy_per_link;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                    }
                }
                #[test]
                fn nondimensional_helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let nondimensional_helmholtz_free_energy = model.legendre.nondimensional_helmholtz_free_energy(&nondimensional_force, temperature);
                        let nondimensional_helmholtz_free_energy_per_link = model.legendre.nondimensional_helmholtz_free_energy_per_link(&nondimensional_force, temperature);
                        let residual_abs = &nondimensional_helmholtz_free_energy/(model.number_of_links as f64) - &nondimensional_helmholtz_free_energy_per_link;
                        let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                    }
                }
                #[test]
                fn nondimensional_relative_helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let nondimensional_relative_helmholtz_free_energy = model.legendre.nondimensional_relative_helmholtz_free_energy(&nondimensional_force);
                        let nondimensional_relative_helmholtz_free_energy_per_link = model.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force);
                        let residual_abs = &nondimensional_relative_helmholtz_free_energy/(model.number_of_links as f64) - &nondimensional_relative_helmholtz_free_energy_per_link;
                        let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                    }
                }
            }
        }
        mod relative
        {
            use super::*;
            use rand::Rng;
            use crate::physics::single_chain::test::Parameters;
            #[test]
            fn gibbs_free_energy()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let gibbs_free_energy = model.gibbs_free_energy(&force, temperature);
                    let gibbs_free_energy_0 = model.gibbs_free_energy(&parameters.zero, temperature);
                    let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&force, temperature);
                    let residual_abs = &gibbs_free_energy - &gibbs_free_energy_0 - &relative_gibbs_free_energy;
                    let residual_rel = &residual_abs/&gibbs_free_energy_0;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn gibbs_free_energy_per_link()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&force, temperature);
                    let gibbs_free_energy_per_link_0 = model.gibbs_free_energy_per_link(&parameters.zero, temperature);
                    let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&force, temperature);
                    let residual_abs = &gibbs_free_energy_per_link - &gibbs_free_energy_per_link_0 - &relative_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&gibbs_free_energy_per_link_0;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn nondimensional_gibbs_free_energy()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_force, temperature);
                    let nondimensional_gibbs_free_energy_0 = model.nondimensional_gibbs_free_energy(&parameters.zero, temperature);
                    let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
                    let residual_abs = &nondimensional_gibbs_free_energy - &nondimensional_gibbs_free_energy_0 - &nondimensional_relative_gibbs_free_energy;
                    let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_0;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn nondimensional_gibbs_free_energy_per_link()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, temperature);
                    let nondimensional_gibbs_free_energy_per_link_0 = model.nondimensional_gibbs_free_energy_per_link(&parameters.zero, temperature);
                    let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
                    let residual_abs = &nondimensional_gibbs_free_energy_per_link - &nondimensional_gibbs_free_energy_per_link_0 - &nondimensional_relative_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link_0;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
            mod legendre
            {
                use super::*;
                use rand::Rng;
                use crate::physics::single_chain::test::Parameters;
                #[test]
                fn helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let helmholtz_free_energy = model.legendre.helmholtz_free_energy(&force, temperature);
                        let helmholtz_free_energy_0 = model.legendre.helmholtz_free_energy(&parameters.zero, temperature);
                        let relative_helmholtz_free_energy = model.legendre.relative_helmholtz_free_energy(&force, temperature);
                        let residual_abs = &helmholtz_free_energy - &helmholtz_free_energy_0 - &relative_helmholtz_free_energy;
                        let residual_rel = &residual_abs/&helmholtz_free_energy_0;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                    }
                }
                #[test]
                fn helmholtz_free_energy_per_link()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let helmholtz_free_energy_per_link = model.legendre.helmholtz_free_energy_per_link(&force, temperature);
                        let helmholtz_free_energy_per_link_0 = model.legendre.helmholtz_free_energy_per_link(&parameters.zero, temperature);
                        let relative_helmholtz_free_energy_per_link = model.legendre.relative_helmholtz_free_energy_per_link(&force, temperature);
                        let residual_abs = &helmholtz_free_energy_per_link - &helmholtz_free_energy_per_link_0 - &relative_helmholtz_free_energy_per_link;
                        let residual_rel = &residual_abs/&helmholtz_free_energy_per_link_0;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                    }
                }
                #[test]
                fn nondimensional_helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let nondimensional_helmholtz_free_energy = model.legendre.nondimensional_helmholtz_free_energy(&nondimensional_force, temperature);
                        let nondimensional_helmholtz_free_energy_0 = model.legendre.nondimensional_helmholtz_free_energy(&parameters.zero, temperature);
                        let nondimensional_relative_helmholtz_free_energy = model.legendre.nondimensional_relative_helmholtz_free_energy(&nondimensional_force);
                        let residual_abs = &nondimensional_helmholtz_free_energy - &nondimensional_helmholtz_free_energy_0 - &nondimensional_relative_helmholtz_free_energy;
                        let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_0;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                    }
                }
                #[test]
                fn nondimensional_helmholtz_free_energy_per_link()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let nondimensional_helmholtz_free_energy_per_link = model.legendre.nondimensional_helmholtz_free_energy_per_link(&nondimensional_force, temperature);
                        let nondimensional_helmholtz_free_energy_per_link_0 = model.legendre.nondimensional_helmholtz_free_energy_per_link(&parameters.zero, temperature);
                        let nondimensional_relative_helmholtz_free_energy_per_link = model.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force);
                        let residual_abs = &nondimensional_helmholtz_free_energy_per_link - &nondimensional_helmholtz_free_energy_per_link_0 - &nondimensional_relative_helmholtz_free_energy_per_link;
                        let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link_0;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                    }
                }
            }
        }
        mod slope
        {
            use super::*;
            use rand::Rng;
            use crate::physics::single_chain::test::Parameters;
            #[test]
            fn end_to_end_length()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let small_force = parameters.small*&force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = model.end_to_end_length(&small_force, temperature);
                    let end_to_end_length_from_initial_slope = &small_force/BOLTZMANN_CONSTANT/temperature*(number_of_links as f64)*link_length.powf(2.0)/3.0;
                    let residual_abs = &end_to_end_length - &end_to_end_length_from_initial_slope;
                    let residual_rel = &residual_abs/&end_to_end_length;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                }
            }
            #[test]
            fn end_to_end_length_per_link()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let small_force = parameters.small*&force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length_per_link = model.end_to_end_length_per_link(&small_force, temperature);
                    let end_to_end_length_per_link_from_initial_slope = &small_force/BOLTZMANN_CONSTANT/temperature*link_length.powf(2.0)/3.0;
                    let residual_abs = &end_to_end_length_per_link - &end_to_end_length_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&end_to_end_length_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                }
            }
            #[test]
            fn nondimensional_end_to_end_length()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_force = parameters.small*&nondimensional_force;
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&small_nondimensional_force);
                    let nondimensional_end_to_end_length_from_initial_slope = &small_nondimensional_force*(number_of_links as f64)/3.0;
                    let residual_abs = &nondimensional_end_to_end_length - &nondimensional_end_to_end_length_from_initial_slope;
                    let residual_rel = &residual_abs/&nondimensional_end_to_end_length;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                }
            }
            #[test]
            fn nondimensional_end_to_end_length_per_link()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_force = parameters.small*&nondimensional_force;
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&small_nondimensional_force);
                    let nondimensional_end_to_end_length_per_link_from_initial_slope = &small_nondimensional_force/3.0;
                    let residual_abs = &nondimensional_end_to_end_length_per_link - &nondimensional_end_to_end_length_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                }
            }
            // #[test]
            // fn gibbs_free_energy()
            // {
            // }
            // #[test]
            // fn gibbs_free_energy_per_link()
            // {
            // }
            #[test]
            fn relative_gibbs_free_energy()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let small_force = parameters.small*&force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&small_force, temperature);
                    let relative_gibbs_free_energy_from_initial_slope = -(small_force/BOLTZMANN_CONSTANT/temperature*link_length).powf(2.0)*(number_of_links as f64)*BOLTZMANN_CONSTANT*temperature/6.0;
                    let residual_abs = &relative_gibbs_free_energy - &relative_gibbs_free_energy_from_initial_slope;
                    let residual_rel = &residual_abs/&relative_gibbs_free_energy;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                }
            }
            #[test]
            fn relative_gibbs_free_energy_per_link()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let small_force = parameters.small*&force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&small_force, temperature);
                    let relative_gibbs_free_energy_per_link_from_initial_slope = -(small_force/BOLTZMANN_CONSTANT/temperature*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature/6.0;
                    let residual_abs = &relative_gibbs_free_energy_per_link - &relative_gibbs_free_energy_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&relative_gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                }
            }
            // #[test]
            // fn nondimensional_gibbs_free_energy()
            // {
            // }
            // #[test]
            // fn nondimensional_gibbs_free_energy_per_link()
            // {
            // }
            #[test]
            fn nondimensional_relative_gibbs_free_energy()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_force = parameters.small*&nondimensional_force;
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&small_nondimensional_force);
                    let nondimensional_relative_gibbs_free_energy_from_initial_slope = -small_nondimensional_force.powf(2.0)*(number_of_links as f64)/6.0;
                    let residual_abs = &nondimensional_relative_gibbs_free_energy - &nondimensional_relative_gibbs_free_energy_from_initial_slope;
                    let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                }
            }
            #[test]
            fn nondimensional_relative_gibbs_free_energy_per_link()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_force = parameters.small*&nondimensional_force;
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&small_nondimensional_force);
                    let nondimensional_relative_gibbs_free_energy_per_link_from_initial_slope = -small_nondimensional_force.powf(2.0)/6.0;
                    let residual_abs = &nondimensional_relative_gibbs_free_energy_per_link - &nondimensional_relative_gibbs_free_energy_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                }
            }
            mod legendre
            {
                use super::*;
                use rand::Rng;
                use crate::physics::single_chain::test::Parameters;
                // #[test]
                // fn helmholtz_free_energy()
                // {
                // }
                // #[test]
                // fn helmholtz_free_energy_per_link()
                // {
                // }
                #[test]
                fn relative_helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                        let small_force = parameters.small*&force;
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let relative_helmholtz_free_energy = model.legendre.relative_helmholtz_free_energy(&small_force, temperature);
                        let relative_helmholtz_free_energy_from_initial_slope = (small_force/BOLTZMANN_CONSTANT/temperature*link_length).powf(2.0)*(number_of_links as f64)*BOLTZMANN_CONSTANT*temperature/6.0;
                        let residual_abs = &relative_helmholtz_free_energy - &relative_helmholtz_free_energy_from_initial_slope;
                        let residual_rel = &residual_abs/&relative_helmholtz_free_energy;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                    }
                }
                #[test]
                fn relative_helmholtz_free_energy_per_link()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                        let small_force = parameters.small*&force;
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let relative_helmholtz_free_energy_per_link = model.legendre.relative_helmholtz_free_energy_per_link(&small_force, temperature);
                        let relative_helmholtz_free_energy_per_link_from_initial_slope = (small_force/BOLTZMANN_CONSTANT/temperature*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature/6.0;
                        let residual_abs = &relative_helmholtz_free_energy_per_link - &relative_helmholtz_free_energy_per_link_from_initial_slope;
                        let residual_rel = &residual_abs/&relative_helmholtz_free_energy_per_link;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                    }
                }
                // #[test]
                // fn nondimensional_helmholtz_free_energy()
                // {
                // }
                // #[test]
                // fn nondimensional_helmholtz_free_energy_per_link()
                // {
                // }
                #[test]
                fn nondimensional_relative_helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                        let small_nondimensional_force = parameters.small*&nondimensional_force;
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let nondimensional_relative_helmholtz_free_energy = model.legendre.nondimensional_relative_helmholtz_free_energy(&small_nondimensional_force);
                        let nondimensional_relative_helmholtz_free_energy_from_initial_slope = small_nondimensional_force.powf(2.0)*(number_of_links as f64)/6.0;
                        let residual_abs = &nondimensional_relative_helmholtz_free_energy - &nondimensional_relative_helmholtz_free_energy_from_initial_slope;
                        let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                    }
                }
                #[test]
                fn nondimensional_relative_helmholtz_free_energy_per_link()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                        let small_nondimensional_force = parameters.small*&nondimensional_force;
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let nondimensional_relative_helmholtz_free_energy_per_link = model.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&small_nondimensional_force);
                        let nondimensional_relative_helmholtz_free_energy_per_link_from_initial_slope = small_nondimensional_force.powf(2.0)/6.0;
                        let residual_abs = &nondimensional_relative_helmholtz_free_energy_per_link - &nondimensional_relative_helmholtz_free_energy_per_link_from_initial_slope;
                        let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
                        assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                        assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                    }
                }
            }
        }
    }
}
pub(crate) use isotensional;