#![cfg(test)]

pub static TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED: f64 = 1e-2;

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
    pub end_to_end_length_reference: f64,
    pub end_to_end_length_scale: f64,
    pub nondimensional_end_to_end_length_per_link_reference: f64,
    pub nondimensional_end_to_end_length_per_link_scale: f64,
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
    pub large_number_of_links: u16,
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
            number_of_loops: 888,
            abs_tol_for_close: 1e-4,
            rel_tol_for_close: 1e-3,
            abs_tol_for_equals: 1e-5,
            rel_tol_for_equals: 1e-4,
            force_reference: 1e0,
            force_scale: 1e0,
            nondimensional_force_reference: 1e0,
            nondimensional_force_scale: 1e0,
            end_to_end_length_reference: 1e-1,
            end_to_end_length_scale: 1e-1,
            nondimensional_end_to_end_length_per_link_reference: 5e-1,
            nondimensional_end_to_end_length_per_link_scale: 49e-2,
            temperature_reference: 3e2,
            temperature_scale: 1e0,
            hinge_mass_reference: 1e0,
            hinge_mass_scale: 1e0,
            link_length_reference: 1e0,
            link_length_scale: 1e0,
            default_link_length: 1e0,
            minimum_number_of_links: 6,
            maximum_number_of_links: 25,
            default_number_of_links: 8,
            large_number_of_links: 25,
            default_hinge_mass: 1e0,
            zero: 1e-8,
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
        // mod init
        // {
        //     use super::*;
        //     use crate::physics::single_chain::test::Parameters;
        //     #[test]
        //     fn number_of_configurations()
        //     {}
        // }
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
        mod legendre
        {
            use super::*;
            use rand::Rng;
            use crate::physics::single_chain::test::Parameters;
            use crate::physics::single_chain::test::TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED;
            #[test]
            fn force()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force_in = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let end_to_end_length = model.isotensional.end_to_end_length(&force_in, temperature);
                    let force_out = model.isometric.legendre.force(&end_to_end_length, temperature);
                    let residual_abs = &force_in - &force_out;
                    let residual_rel = &residual_abs/&force_in;
                    assert!(residual_abs.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
                    assert!(residual_rel.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
                }
            }
            #[test]
            fn nondimensional_force()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force_in = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_end_to_end_length_per_link = model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force_in);
                    let nondimensional_force_out = model.isometric.legendre.nondimensional_force(&nondimensional_end_to_end_length_per_link);
                    let residual_abs = &nondimensional_force_in - &nondimensional_force_out;
                    let residual_rel = &residual_abs/&nondimensional_force_in;
                    assert!(residual_abs.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
                    assert!(residual_rel.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
                }
            }
            // #[test]
            // fn end_to_end_length()
            // {}
            // #[test]
            // fn end_to_end_length_per_link()
            // {}
            // #[test]
            // fn nondimensional_end_to_end_length()
            // {}
            // #[test]
            // fn nondimensional_end_to_end_length_per_link()
            // {}
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy_isotensional_legendre = model.isotensional.legendre.helmholtz_free_energy(&force, temperature);
                    let end_to_end_length = model.isotensional.end_to_end_length(&force, temperature);
                    let helmholtz_free_energy_isometric_legendre = model.isometric.legendre.helmholtz_free_energy(&end_to_end_length, temperature);
                    let residual_abs = &helmholtz_free_energy_isotensional_legendre - &helmholtz_free_energy_isometric_legendre;
                    let residual_rel = &residual_abs/&helmholtz_free_energy_isotensional_legendre;
                    assert!(residual_abs.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
                    assert!(residual_rel.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy_per_link_isotensional_legendre = model.isotensional.legendre.helmholtz_free_energy_per_link(&force, temperature);
                    let end_to_end_length = model.isotensional.end_to_end_length(&force, temperature);
                    let helmholtz_free_energy_per_link_isometric_legendre = model.isometric.legendre.helmholtz_free_energy_per_link(&end_to_end_length, temperature);
                    let residual_abs = &helmholtz_free_energy_per_link_isotensional_legendre - &helmholtz_free_energy_per_link_isometric_legendre;
                    let residual_rel = &residual_abs/&helmholtz_free_energy_per_link_isotensional_legendre;
                    assert!(residual_abs.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
                    assert!(residual_rel.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_helmholtz_free_energy_isotensional_legendre = model.isotensional.legendre.relative_helmholtz_free_energy(&force, temperature);
                    let end_to_end_length = model.isotensional.end_to_end_length(&force, temperature);
                    let relative_helmholtz_free_energy_isometric_legendre = model.isometric.legendre.relative_helmholtz_free_energy(&end_to_end_length, temperature);
                    let residual_abs = &relative_helmholtz_free_energy_isotensional_legendre - &relative_helmholtz_free_energy_isometric_legendre;
                    let residual_rel = &residual_abs/&relative_helmholtz_free_energy_isotensional_legendre;
                    assert!(residual_abs.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
                    assert!(residual_rel.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_helmholtz_free_energy_per_link_isotensional_legendre = model.isotensional.legendre.relative_helmholtz_free_energy_per_link(&force, temperature);
                    let end_to_end_length = model.isotensional.end_to_end_length(&force, temperature);
                    let relative_helmholtz_free_energy_per_link_isometric_legendre = model.isometric.legendre.relative_helmholtz_free_energy_per_link(&end_to_end_length, temperature);
                    let residual_abs = &relative_helmholtz_free_energy_per_link_isotensional_legendre - &relative_helmholtz_free_energy_per_link_isometric_legendre;
                    let residual_rel = &residual_abs/&relative_helmholtz_free_energy_per_link_isotensional_legendre;
                    assert!(residual_abs.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
                    assert!(residual_rel.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_helmholtz_free_energy_isotensional_legendre = model.isotensional.legendre.nondimensional_helmholtz_free_energy(&nondimensional_force, temperature);
                    let nondimensional_end_to_end_length_per_link = model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                    let nondimensional_helmholtz_free_energy_isometric_legendre = model.isometric.legendre.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, temperature);
                    let residual_abs = &nondimensional_helmholtz_free_energy_isotensional_legendre - &nondimensional_helmholtz_free_energy_isometric_legendre;
                    let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_isotensional_legendre;
                    assert!(residual_abs.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
                    assert!(residual_rel.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_helmholtz_free_energy_per_link_isotensional_legendre = model.isotensional.legendre.nondimensional_helmholtz_free_energy_per_link(&nondimensional_force, temperature);
                    let nondimensional_end_to_end_length_per_link = model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                    let nondimensional_helmholtz_free_energy_per_link_isometric_legendre = model.isometric.legendre.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, temperature);
                    let residual_abs = &nondimensional_helmholtz_free_energy_per_link_isotensional_legendre - &nondimensional_helmholtz_free_energy_per_link_isometric_legendre;
                    let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link_isotensional_legendre;
                    assert!(residual_abs.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
                    assert!(residual_rel.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_relative_helmholtz_free_energy_isotensional_legendre = model.isotensional.legendre.nondimensional_relative_helmholtz_free_energy(&nondimensional_force);
                    let nondimensional_end_to_end_length_per_link = model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                    let nondimensional_relative_helmholtz_free_energy_isometric_legendre = model.isometric.legendre.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link);
                    let residual_abs = &nondimensional_relative_helmholtz_free_energy_isotensional_legendre - &nondimensional_relative_helmholtz_free_energy_isometric_legendre;
                    let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_isotensional_legendre;
                    assert!(residual_abs.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
                    assert!(residual_rel.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_relative_helmholtz_free_energy_per_link_isotensional_legendre = model.isotensional.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force);
                    let nondimensional_end_to_end_length_per_link = model.isotensional.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                    let nondimensional_relative_helmholtz_free_energy_per_link_isometric_legendre = model.isometric.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
                    let residual_abs = &nondimensional_relative_helmholtz_free_energy_per_link_isotensional_legendre - &nondimensional_relative_helmholtz_free_energy_per_link_isometric_legendre;
                    let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link_isotensional_legendre;
                    assert!(residual_abs.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
                    assert!(residual_rel.abs() <= TEMPORARY_REDUCED_TOL_FOR_INV_LANG_RELATED);
                }
            }
            // #[test]
            // fn gibbs_free_energy()
            // {}
            // #[test]
            // fn gibbs_free_energy_per_link()
            // {}
            // #[test]
            // fn relative_gibbs_free_energy()
            // {}
            // #[test]
            // fn relative_gibbs_free_energy_per_link()
            // {}
            // #[test]
            // fn nondimensional_gibbs_free_energy()
            // {}
            // #[test]
            // fn nondimensional_gibbs_free_energy_per_link()
            // {}
            // #[test]
            // fn nondimensional_relative_gibbs_free_energy()
            // {}
            // #[test]
            // fn nondimensional_relative_gibbs_free_energy_per_link()
            // {}
            mod thermodynamic_limit
            {
                use super::*;
                use rand::Rng;
                use crate::physics::BOLTZMANN_CONSTANT;
                use crate::physics::single_chain::test::Parameters;
                #[test]
                fn force()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = parameters.large_number_of_links;
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let force = model.isometric.force(&end_to_end_length, temperature);
                        let force_legendre = model.isometric.legendre.force(&end_to_end_length, temperature);
                        let residual_abs = &force_legendre - &force;
                        let residual_rel = &residual_abs/&force;
                        assert!(residual_abs.abs() <= 3e-1*(number_of_links as f64).sqrt());
                        assert!(residual_rel.abs() <= 3e-1);
                    }
                }
                #[test]
                fn nondimensional_force()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = parameters.large_number_of_links;
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                        let nondimensional_force = model.isometric.nondimensional_force(&nondimensional_end_to_end_length_per_link);
                        let nondimensional_force_legendre = model.isometric.legendre.nondimensional_force(&nondimensional_end_to_end_length_per_link);
                        let residual_abs = &nondimensional_force_legendre - &nondimensional_force;
                        let residual_rel = &residual_abs/&nondimensional_force;
                        assert!(residual_abs.abs() <= 3e-1*(number_of_links as f64).sqrt());
                        assert!(residual_rel.abs() <= 3e-1);
                    }
                }
                // #[test]
                // fn end_to_end_length()
                // {}
                // #[test]
                // fn end_to_end_length_per_link()
                // {}
                // #[test]
                // fn nondimensional_end_to_end_length()
                // {}
                // #[test]
                // fn nondimensional_end_to_end_length_per_link()
                // {}
                #[test]
                fn helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = parameters.large_number_of_links;
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let helmholtz_free_energy = model.isometric.helmholtz_free_energy(&end_to_end_length, temperature);
                        let helmholtz_free_energy_legendre = model.isometric.legendre.helmholtz_free_energy(&end_to_end_length, temperature);
                        let residual_abs = &helmholtz_free_energy_legendre - &helmholtz_free_energy;
                        let residual_rel = &residual_abs/&helmholtz_free_energy;
                        // assert!(residual_abs.abs() <= 3e-1*BOLTZMANN_CONSTANT*temperature*(number_of_links as f64).sqrt());
                        assert!(residual_rel.abs() <= 3e-1);
                    }
                }
                #[test]
                fn helmholtz_free_energy_per_link()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = parameters.large_number_of_links;
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let helmholtz_free_energy_per_link = model.isometric.helmholtz_free_energy_per_link(&end_to_end_length, temperature);
                        let helmholtz_free_energy_per_link_legendre = model.isometric.legendre.helmholtz_free_energy_per_link(&end_to_end_length, temperature);
                        let residual_abs = &helmholtz_free_energy_per_link_legendre - &helmholtz_free_energy_per_link;
                        let residual_rel = &residual_abs/&helmholtz_free_energy_per_link;
                        // assert!(residual_abs.abs() <= 3e-1*BOLTZMANN_CONSTANT*temperature/(number_of_links as f64).sqrt());
                        assert!(residual_rel.abs() <= 3e-1);
                    }
                }
                #[test]
                fn relative_helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = parameters.large_number_of_links;
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let relative_helmholtz_free_energy = model.isometric.relative_helmholtz_free_energy(&end_to_end_length, temperature);
                        let relative_helmholtz_free_energy_legendre = model.isometric.legendre.relative_helmholtz_free_energy(&end_to_end_length, temperature);
                        let residual_abs = &relative_helmholtz_free_energy_legendre - &relative_helmholtz_free_energy;
                        // let residual_rel = &residual_abs/&relative_helmholtz_free_energy;
                        assert!(residual_abs.abs() <= 3e-1*BOLTZMANN_CONSTANT*temperature/(number_of_links as f64).sqrt());
                        // assert!(residual_rel.abs() <= 3e-1);
                    }
                }
                #[test]
                fn relative_helmholtz_free_energy_per_link()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = parameters.large_number_of_links;
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let relative_helmholtz_free_energy_per_link = model.isometric.relative_helmholtz_free_energy_per_link(&end_to_end_length, temperature);
                        let relative_helmholtz_free_energy_per_link_legendre = model.isometric.legendre.relative_helmholtz_free_energy_per_link(&end_to_end_length, temperature);
                        let residual_abs = &relative_helmholtz_free_energy_per_link_legendre - &relative_helmholtz_free_energy_per_link;
                        // let residual_rel = &residual_abs/&relative_helmholtz_free_energy_per_link;
                        assert!(residual_abs.abs() <= 3e-1*BOLTZMANN_CONSTANT*temperature/(number_of_links as f64).sqrt());
                        // assert!(residual_rel.abs() <= 3e-1);
                    }
                }
                #[test]
                fn nondimensional_helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = parameters.large_number_of_links;
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let nondimensional_helmholtz_free_energy = model.isometric.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, temperature);
                        let nondimensional_helmholtz_free_energy_legendre = model.isometric.legendre.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, temperature);
                        let residual_abs = &nondimensional_helmholtz_free_energy_legendre - &nondimensional_helmholtz_free_energy;
                        let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy;
                        // assert!(residual_abs.abs() <= 3e-1*(number_of_links as f64).sqrt());
                        assert!(residual_rel.abs() <= 3e-1);
                    }
                }
                #[test]
                fn nondimensional_helmholtz_free_energy_per_link()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = parameters.large_number_of_links;
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                        let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                        let nondimensional_helmholtz_free_energy_per_link = model.isometric.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, temperature);
                        let nondimensional_helmholtz_free_energy_per_link_legendre = model.isometric.legendre.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, temperature);
                        let residual_abs = &nondimensional_helmholtz_free_energy_per_link_legendre - &nondimensional_helmholtz_free_energy_per_link;
                        let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link;
                        // assert!(residual_abs.abs() <= 3e-1/(number_of_links as f64).sqrt());
                        assert!(residual_rel.abs() <= 3e-1);
                    }
                }
                #[test]
                fn nondimensional_relative_helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = parameters.large_number_of_links;
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                        let nondimensional_relative_helmholtz_free_energy = model.isometric.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link);
                        let nondimensional_relative_helmholtz_free_energy_legendre = model.isometric.legendre.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link);
                        let residual_abs = &nondimensional_relative_helmholtz_free_energy_legendre - &nondimensional_relative_helmholtz_free_energy;
                        let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy;
                        assert!(residual_abs.abs() <= 3e-1*(number_of_links as f64).sqrt());
                        assert!(residual_rel.abs() <= 3e-1);
                    }
                }
                #[test]
                fn nondimensional_relative_helmholtz_free_energy_per_link()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = parameters.large_number_of_links;
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                        let nondimensional_relative_helmholtz_free_energy_per_link = model.isometric.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
                        let nondimensional_relative_helmholtz_free_energy_per_link_legendre = model.isometric.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
                        let residual_abs = &nondimensional_relative_helmholtz_free_energy_per_link_legendre - &nondimensional_relative_helmholtz_free_energy_per_link;
                        let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
                        assert!(residual_abs.abs() <= 3e-1/(number_of_links as f64).sqrt());
                        assert!(residual_rel.abs() <= 3e-1);
                    }
                }
                #[test]
                fn equilibrium_distribution()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                        let equilibrium_distribution = model.isometric.equilibrium_distribution(&end_to_end_length);
                        let equilibrium_distribution_legendre = model.isometric.legendre.equilibrium_distribution(&end_to_end_length);
                        let residual_abs = &equilibrium_distribution_legendre - &equilibrium_distribution;
                        let residual_rel = &residual_abs/&equilibrium_distribution;
                        assert!(residual_abs.abs() <= 3e-1/(number_of_links as f64).sqrt() || residual_rel.abs() <= 3e-1);
                    }
                }
                #[test]
                fn nondimensional_equilibrium_distribution()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                        let nondimensional_equilibrium_distribution = model.isometric.nondimensional_equilibrium_distribution(&nondimensional_end_to_end_length_per_link);
                        let nondimensional_equilibrium_distribution_legendre = model.isometric.legendre.nondimensional_equilibrium_distribution(&nondimensional_end_to_end_length_per_link);
                        let residual_abs = &nondimensional_equilibrium_distribution_legendre - &nondimensional_equilibrium_distribution;
                        let residual_rel = &residual_abs/&nondimensional_equilibrium_distribution;
                        assert!(residual_abs.abs() <= 3e-1/(number_of_links as f64).sqrt() || residual_rel.abs() <= 3e-1);
                    }
                }
                #[test]
                fn equilibrium_radial_distribution()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                        let equilibrium_radial_distribution = model.isometric.equilibrium_radial_distribution(&end_to_end_length);
                        let equilibrium_radial_distribution_legendre = model.isometric.legendre.equilibrium_radial_distribution(&end_to_end_length);
                        let residual_abs = &equilibrium_radial_distribution_legendre - &equilibrium_radial_distribution;
                        let residual_rel = &residual_abs/&equilibrium_radial_distribution;
                        assert!(residual_abs.abs() <= 3e-1/(number_of_links as f64).sqrt() || residual_rel.abs() <= 3e-1);
                    }
                }
                #[test]
                fn nondimensional_equilibrium_radial_distribution()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                        let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                        let model = <$model>::init(number_of_links, link_length, hinge_mass);
                        let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                        let nondimensional_equilibrium_radial_distribution = model.isometric.nondimensional_equilibrium_radial_distribution(&nondimensional_end_to_end_length_per_link);
                        let nondimensional_equilibrium_radial_distribution_legendre = model.isometric.legendre.nondimensional_equilibrium_radial_distribution(&nondimensional_end_to_end_length_per_link);
                        let residual_abs = &nondimensional_equilibrium_radial_distribution_legendre - &nondimensional_equilibrium_radial_distribution;
                        let residual_rel = &residual_abs/&nondimensional_equilibrium_radial_distribution;
                        // assert!(residual_abs.abs() <= 3e-1/(number_of_links as f64).sqrt() || residual_rel.abs() <= 3e-1);
                        assert!(residual_abs.abs() <= 3e1/(number_of_links as f64).sqrt() || residual_rel.abs() <= 3e1);
                    }
                }
                // #[test]
                // fn gibbs_free_energy()
                // {}
                // #[test]
                // fn gibbs_free_energy_per_link()
                // {}
                // #[test]
                // fn relative_gibbs_free_energy()
                // {}
                // #[test]
                // fn relative_gibbs_free_energy_per_link()
                // {}
                // #[test]
                // fn nondimensional_gibbs_free_energy()
                // {}
                // #[test]
                // fn nondimensional_gibbs_free_energy_per_link()
                // {}
                // #[test]
                // fn nondimensional_relative_gibbs_free_energy()
                // {}
                // #[test]
                // fn nondimensional_relative_gibbs_free_energy_per_link()
                // {}
            }
        }
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
            fn init()
            {
                let parameters = Parameters::default();
                let _ = <$model>::init(parameters.default_number_of_links, parameters.default_link_length, parameters.default_hinge_mass);
            }
        }
        mod normalization
        {
            use super::*;
            use rand::Rng;
            use crate::math::integrate;
            use crate::physics::single_chain::test::Parameters;
            #[test]
            fn equilibrium_distribution()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let integrand = |end_to_end_length: f64| 4.0*PI*end_to_end_length.powf(2.0)*model.equilibrium_distribution(&end_to_end_length);
                    let integral = integrate(integrand, 0.0, model.contour_length, 100);
                    assert!((integral - 1.0).abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn nondimensional_equilibrium_distribution()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let integrand = |nondimensional_end_to_end_length_per_link: f64| 4.0*PI*nondimensional_end_to_end_length_per_link.powf(2.0)*model.nondimensional_equilibrium_distribution(&nondimensional_end_to_end_length_per_link);
                    let integral = integrate(integrand, 0.0, 1.0, 100);
                    assert!((integral - 1.0).abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn equilibrium_radial_distribution()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let integrand = |end_to_end_length: f64| model.equilibrium_radial_distribution(&end_to_end_length);
                    let integral = integrate(integrand, 0.0, model.contour_length, 100);
                    assert!((integral - 1.0).abs() <= parameters.rel_tol_for_equals);
                }
            }
            #[test]
            fn nondimensional_equilibrium_radial_distribution()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let integrand = |nondimensional_end_to_end_length_per_link: f64| model.nondimensional_equilibrium_radial_distribution(&nondimensional_end_to_end_length_per_link);
                    let integral = integrate(integrand, 0.0, 1.0, 100);
                    assert!((integral - 1.0).abs() <= parameters.rel_tol_for_equals);
                }
            }
        }
        mod nondimensional
        {
            use super::*;
            use rand::Rng;
            use crate::physics::single_chain::test::Parameters;
            #[test]
            fn force()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let force = model.force(&end_to_end_length, temperature);
                    let nondimensional_end_to_end_length_per_link = end_to_end_length/(number_of_links as f64)/link_length;
                    let nondimensional_force = model.nondimensional_force(&nondimensional_end_to_end_length_per_link);
                    let residual_abs = &force/BOLTZMANN_CONSTANT/temperature*link_length - &nondimensional_force;
                    let residual_rel = &residual_abs/&nondimensional_force;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy = model.helmholtz_free_energy(&end_to_end_length, temperature);
                    let nondimensional_end_to_end_length_per_link = end_to_end_length/(number_of_links as f64)/link_length;
                    let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, temperature);
                    let residual_abs = &helmholtz_free_energy/BOLTZMANN_CONSTANT/temperature - &nondimensional_helmholtz_free_energy;
                    let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy;
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&end_to_end_length, temperature);
                    let nondimensional_end_to_end_length_per_link = end_to_end_length/(number_of_links as f64)/link_length;
                    let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, temperature);
                    let residual_abs = &helmholtz_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_helmholtz_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link;
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&end_to_end_length, temperature);
                    let nondimensional_end_to_end_length_per_link = end_to_end_length/(number_of_links as f64)/link_length;
                    let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link);
                    let residual_abs = &relative_helmholtz_free_energy/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_helmholtz_free_energy;
                    let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy;
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&end_to_end_length, temperature);
                    let nondimensional_end_to_end_length_per_link = end_to_end_length/(number_of_links as f64)/link_length;
                    let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
                    let residual_abs = &relative_helmholtz_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_helmholtz_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
        }
        mod per_link
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy = model.helmholtz_free_energy(&end_to_end_length, temperature);
                    let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&end_to_end_length, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&end_to_end_length, temperature);
                    let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&end_to_end_length, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, temperature);
                    let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link);
                    let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
                    let residual_abs = &nondimensional_relative_helmholtz_free_energy/(model.number_of_links as f64) - &nondimensional_relative_helmholtz_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
        }
        mod relative
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy = model.helmholtz_free_energy(&end_to_end_length, temperature);
                    let helmholtz_free_energy_0 = model.helmholtz_free_energy(&parameters.zero, temperature);
                    let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&end_to_end_length, temperature);
                    let residual_abs = &helmholtz_free_energy - &helmholtz_free_energy_0 - &relative_helmholtz_free_energy;
                    let residual_rel = &residual_abs/&helmholtz_free_energy_0;
                    // assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&end_to_end_length, temperature);
                    let helmholtz_free_energy_per_link_0 = model.helmholtz_free_energy_per_link(&parameters.zero, temperature);
                    let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&end_to_end_length, temperature);
                    let residual_abs = &helmholtz_free_energy_per_link - &helmholtz_free_energy_per_link_0 - &relative_helmholtz_free_energy_per_link;
                    let residual_rel = &residual_abs/&helmholtz_free_energy_per_link_0;
                    // assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link, temperature);
                    let nondimensional_helmholtz_free_energy_0 = model.nondimensional_helmholtz_free_energy(&parameters.zero, temperature);
                    let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_end_to_end_length_per_link);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link, temperature);
                    let nondimensional_helmholtz_free_energy_per_link_0 = model.nondimensional_helmholtz_free_energy_per_link(&parameters.zero, temperature);
                    let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_end_to_end_length_per_link);
                    let residual_abs = &nondimensional_helmholtz_free_energy_per_link - &nondimensional_helmholtz_free_energy_per_link_0 - &nondimensional_relative_helmholtz_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link_0;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
        }
        // mod slope
        // {if do, move helmholtz from isometricLegendre here, and put in gibbs in isometricLegendre}
    }
}
pub(crate) use isometric;

macro_rules! isometricLegendre
{
    ( $model:ty ) =>
    {
        mod nondimensional_legendre
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let gibbs_free_energy = model.gibbs_free_energy(&end_to_end_length, temperature);
                    let nondimensional_end_to_end_length_per_link = end_to_end_length/(number_of_links as f64)/link_length;
                    let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_end_to_end_length_per_link, temperature);
                    let residual_abs = &gibbs_free_energy/BOLTZMANN_CONSTANT/temperature - &nondimensional_gibbs_free_energy;
                    let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy;
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&end_to_end_length, temperature);
                    let nondimensional_end_to_end_length_per_link = end_to_end_length/(number_of_links as f64)/link_length;
                    let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_end_to_end_length_per_link, temperature);
                    let residual_abs = &gibbs_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link;
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&end_to_end_length, temperature);
                    let nondimensional_end_to_end_length_per_link = end_to_end_length/(number_of_links as f64)/link_length;
                    let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_end_to_end_length_per_link, temperature);
                    let residual_abs = &relative_gibbs_free_energy/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_gibbs_free_energy;
                    let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy;
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&end_to_end_length, temperature);
                    let nondimensional_end_to_end_length_per_link = end_to_end_length/(number_of_links as f64)/link_length;
                    let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_end_to_end_length_per_link, temperature);
                    let residual_abs = &relative_gibbs_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
        }
        mod per_link_legendre
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let gibbs_free_energy = model.gibbs_free_energy(&end_to_end_length, temperature);
                    let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&end_to_end_length, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&end_to_end_length, temperature);
                    let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&end_to_end_length, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_end_to_end_length_per_link, temperature);
                    let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_end_to_end_length_per_link, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_end_to_end_length_per_link, temperature);
                    let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_end_to_end_length_per_link, temperature);
                    let residual_abs = &nondimensional_relative_gibbs_free_energy/(model.number_of_links as f64) - &nondimensional_relative_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
        }
        mod relative_legendre
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let gibbs_free_energy = model.gibbs_free_energy(&end_to_end_length, temperature);
                    let gibbs_free_energy_0 = model.gibbs_free_energy(&parameters.zero, temperature);
                    let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&end_to_end_length, temperature);
                    let residual_abs = &gibbs_free_energy - &gibbs_free_energy_0 - &relative_gibbs_free_energy;
                    let residual_rel = &residual_abs/&gibbs_free_energy_0;
                    // assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&end_to_end_length, temperature);
                    let gibbs_free_energy_per_link_0 = model.gibbs_free_energy_per_link(&parameters.zero, temperature);
                    let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&end_to_end_length, temperature);
                    let residual_abs = &gibbs_free_energy_per_link - &gibbs_free_energy_per_link_0 - &relative_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&gibbs_free_energy_per_link_0;
                    // assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&nondimensional_end_to_end_length_per_link, temperature);
                    let nondimensional_gibbs_free_energy_0 = model.nondimensional_gibbs_free_energy(&parameters.zero, temperature);
                    let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_end_to_end_length_per_link, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_end_to_end_length_per_link, temperature);
                    let nondimensional_gibbs_free_energy_per_link_0 = model.nondimensional_gibbs_free_energy_per_link(&parameters.zero, temperature);
                    let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_end_to_end_length_per_link, temperature);
                    let residual_abs = &nondimensional_gibbs_free_energy_per_link - &nondimensional_gibbs_free_energy_per_link_0 - &nondimensional_relative_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link_0;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
        }
        mod slope_legendre
        {
            use super::*;
            use rand::Rng;
            use crate::physics::single_chain::test::Parameters;
            #[test]
            fn force()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let small_end_to_end_length = parameters.small*&end_to_end_length;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let force = model.force(&small_end_to_end_length, temperature);
                    let force_from_initial_slope = (small_end_to_end_length/(number_of_links as f64)/link_length)*BOLTZMANN_CONSTANT*temperature*3.0/link_length;
                    let residual_abs = &force - &force_from_initial_slope;
                    let residual_rel = &residual_abs/&force;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                }
            }
            #[test]
            fn nondimensional_force()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
                    let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_end_to_end_length_per_link = parameters.small*&nondimensional_end_to_end_length_per_link;
                    let nondimensional_force = model.nondimensional_force(&small_nondimensional_end_to_end_length_per_link);
                    let nondimensional_force_from_initial_slope = 3.0*small_nondimensional_end_to_end_length_per_link;
                    let residual_abs = &nondimensional_force - &nondimensional_force_from_initial_slope;
                    let residual_rel = &residual_abs/&nondimensional_force;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                }
            }
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let small_end_to_end_length = parameters.small*&end_to_end_length;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy = model.helmholtz_free_energy(&small_end_to_end_length, temperature);
                    let helmholtz_free_energy_from_initial_slope = model.helmholtz_free_energy(&parameters.zero, temperature) + (small_end_to_end_length/(number_of_links as f64)/link_length).powf(2.0)*(number_of_links as f64)*BOLTZMANN_CONSTANT*temperature*1.5;
                    let residual_abs = &helmholtz_free_energy - &helmholtz_free_energy_from_initial_slope;
                    let residual_rel = &residual_abs/&helmholtz_free_energy;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let small_end_to_end_length = parameters.small*&end_to_end_length;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&small_end_to_end_length, temperature);
                    let helmholtz_free_energy_per_link_from_initial_slope = model.helmholtz_free_energy_per_link(&parameters.zero, temperature) + (small_end_to_end_length/(number_of_links as f64)/link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature*1.5;
                    let residual_abs = &helmholtz_free_energy_per_link - &helmholtz_free_energy_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&helmholtz_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let small_end_to_end_length = parameters.small*&end_to_end_length;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&small_end_to_end_length, temperature);
                    let relative_helmholtz_free_energy_from_initial_slope = (small_end_to_end_length/(number_of_links as f64)/link_length).powf(2.0)*(number_of_links as f64)*BOLTZMANN_CONSTANT*temperature*1.5;
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let end_to_end_length = parameters.end_to_end_length_reference + parameters.end_to_end_length_scale*(0.5 - rng.gen::<f64>());
                    let small_end_to_end_length = parameters.small*&end_to_end_length;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&small_end_to_end_length, temperature);
                    let relative_helmholtz_free_energy_per_link_from_initial_slope = (small_end_to_end_length/(number_of_links as f64)/link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature*1.5;
                    let residual_abs = &relative_helmholtz_free_energy_per_link - &relative_helmholtz_free_energy_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&relative_helmholtz_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_end_to_end_length_per_link = parameters.small*&nondimensional_end_to_end_length_per_link;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&small_nondimensional_end_to_end_length_per_link, temperature);
                    let nondimensional_helmholtz_free_energy_from_initial_slope = model.nondimensional_helmholtz_free_energy(&parameters.zero, temperature) + small_nondimensional_end_to_end_length_per_link.powf(2.0)*(number_of_links as f64)*1.5;
                    let residual_abs = &nondimensional_helmholtz_free_energy - &nondimensional_helmholtz_free_energy_from_initial_slope;
                    let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_end_to_end_length_per_link = parameters.small*&nondimensional_end_to_end_length_per_link;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&small_nondimensional_end_to_end_length_per_link, temperature);
                    let nondimensional_helmholtz_free_energy_per_link_from_initial_slope = model.nondimensional_helmholtz_free_energy_per_link(&parameters.zero, temperature) + small_nondimensional_end_to_end_length_per_link.powf(2.0)*1.5;
                    let residual_abs = &nondimensional_helmholtz_free_energy_per_link - &nondimensional_helmholtz_free_energy_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_end_to_end_length_per_link = parameters.small*&nondimensional_end_to_end_length_per_link;
                    let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&small_nondimensional_end_to_end_length_per_link);
                    let nondimensional_relative_helmholtz_free_energy_from_initial_slope = small_nondimensional_end_to_end_length_per_link.powf(2.0)*(number_of_links as f64)*1.5;
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_end_to_end_length_per_link = parameters.nondimensional_end_to_end_length_per_link_reference + parameters.nondimensional_end_to_end_length_per_link_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_end_to_end_length_per_link = parameters.small*&nondimensional_end_to_end_length_per_link;
                    let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&small_nondimensional_end_to_end_length_per_link);
                    let nondimensional_relative_helmholtz_free_energy_per_link_from_initial_slope = small_nondimensional_end_to_end_length_per_link.powf(2.0)*1.5;
                    let residual_abs = &nondimensional_relative_helmholtz_free_energy_per_link - &nondimensional_relative_helmholtz_free_energy_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                }
            }
        }
    }
}
pub(crate) use isometricLegendre;

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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&force, temperature);
                    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                    let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
                    let residual_abs = &relative_gibbs_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
                    let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
                    let residual_abs = &nondimensional_relative_gibbs_free_energy/(model.number_of_links as f64) - &nondimensional_relative_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&nondimensional_force, temperature);
                    let nondimensional_gibbs_free_energy_per_link_0 = model.nondimensional_gibbs_free_energy_per_link(&parameters.zero, temperature);
                    let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
                    let residual_abs = &nondimensional_gibbs_free_energy_per_link - &nondimensional_gibbs_free_energy_per_link_0 - &nondimensional_relative_gibbs_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link_0;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let small_force = parameters.small*&force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let small_force = parameters.small*&force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_force = parameters.small*&nondimensional_force;
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_force = parameters.small*&nondimensional_force;
                    let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&small_nondimensional_force);
                    let nondimensional_end_to_end_length_per_link_from_initial_slope = &small_nondimensional_force/3.0;
                    let residual_abs = &nondimensional_end_to_end_length_per_link - &nondimensional_end_to_end_length_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let small_force = parameters.small*&force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let gibbs_free_energy = model.gibbs_free_energy(&small_force, temperature);
                    let gibbs_free_energy_from_initial_slope = model.gibbs_free_energy(&parameters.zero, temperature) - (small_force/BOLTZMANN_CONSTANT/temperature*link_length).powf(2.0)*(number_of_links as f64)*BOLTZMANN_CONSTANT*temperature/6.0;
                    let residual_abs = &gibbs_free_energy - &gibbs_free_energy_from_initial_slope;
                    let residual_rel = &residual_abs/&gibbs_free_energy;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let small_force = parameters.small*&force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let gibbs_free_energy_per_link = model.gibbs_free_energy_per_link(&small_force, temperature);
                    let gibbs_free_energy_per_link_from_initial_slope = model.gibbs_free_energy_per_link(&parameters.zero, temperature) - (small_force/BOLTZMANN_CONSTANT/temperature*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature/6.0;
                    let residual_abs = &gibbs_free_energy_per_link - &gibbs_free_energy_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let small_force = parameters.small*&force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let small_force = parameters.small*&force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&small_force, temperature);
                    let relative_gibbs_free_energy_per_link_from_initial_slope = -(small_force/BOLTZMANN_CONSTANT/temperature*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature/6.0;
                    let residual_abs = &relative_gibbs_free_energy_per_link - &relative_gibbs_free_energy_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&relative_gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_force = parameters.small*&nondimensional_force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_gibbs_free_energy = model.nondimensional_gibbs_free_energy(&small_nondimensional_force, temperature);
                    let nondimensional_gibbs_free_energy_from_initial_slope = model.nondimensional_gibbs_free_energy(&parameters.zero, temperature) - small_nondimensional_force.powf(2.0)*(number_of_links as f64)/6.0;
                    let residual_abs = &nondimensional_gibbs_free_energy - &nondimensional_gibbs_free_energy_from_initial_slope;
                    let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_force = parameters.small*&nondimensional_force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_gibbs_free_energy_per_link = model.nondimensional_gibbs_free_energy_per_link(&small_nondimensional_force, temperature);
                    let nondimensional_gibbs_free_energy_per_link_from_initial_slope = model.nondimensional_gibbs_free_energy_per_link(&parameters.zero, temperature) - small_nondimensional_force.powf(2.0)/6.0;
                    let residual_abs = &nondimensional_gibbs_free_energy_per_link - &nondimensional_gibbs_free_energy_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&nondimensional_gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_force = parameters.small*&nondimensional_force;
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_force = parameters.small*&nondimensional_force;
                    let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&small_nondimensional_force);
                    let nondimensional_relative_gibbs_free_energy_per_link_from_initial_slope = -small_nondimensional_force.powf(2.0)/6.0;
                    let residual_abs = &nondimensional_relative_gibbs_free_energy_per_link - &nondimensional_relative_gibbs_free_energy_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
                }
            }
        }
    }
}
pub(crate) use isotensional;

macro_rules! isotensionalLegendre
{
    ( $model:ty ) =>
    {
        mod nondimensional_legendre
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy = model.helmholtz_free_energy(&force, temperature);
                    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                    let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_force, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&force, temperature);
                    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                    let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_force);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&force, temperature);
                    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                    let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_force, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&force, temperature);
                    let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                    let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force);
                    let residual_abs = &relative_helmholtz_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_helmholtz_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
        }
        mod per_link_legendre
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy = model.helmholtz_free_energy(&force, temperature);
                    let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&force, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&force, temperature);
                    let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&force, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_force, temperature);
                    let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_force, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_force);
                    let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force);
                    let residual_abs = &nondimensional_relative_helmholtz_free_energy/(model.number_of_links as f64) - &nondimensional_relative_helmholtz_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
        }
        mod relative_legendre
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy = model.helmholtz_free_energy(&force, temperature);
                    let helmholtz_free_energy_0 = model.helmholtz_free_energy(&parameters.zero, temperature);
                    let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&force, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&force, temperature);
                    let helmholtz_free_energy_per_link_0 = model.helmholtz_free_energy_per_link(&parameters.zero, temperature);
                    let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&force, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&nondimensional_force, temperature);
                    let nondimensional_helmholtz_free_energy_0 = model.nondimensional_helmholtz_free_energy(&parameters.zero, temperature);
                    let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&nondimensional_force);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&nondimensional_force, temperature);
                    let nondimensional_helmholtz_free_energy_per_link_0 = model.nondimensional_helmholtz_free_energy_per_link(&parameters.zero, temperature);
                    let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force);
                    let residual_abs = &nondimensional_helmholtz_free_energy_per_link - &nondimensional_helmholtz_free_energy_per_link_0 - &nondimensional_relative_helmholtz_free_energy_per_link;
                    let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link_0;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_equals);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_equals);
                }
            }
        }
        mod slope_legendre
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let small_force = parameters.small*&force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy = model.helmholtz_free_energy(&small_force, temperature);
                    let helmholtz_free_energy_from_initial_slope = model.helmholtz_free_energy(&parameters.zero, temperature) + (small_force/BOLTZMANN_CONSTANT/temperature*link_length).powf(2.0)*(number_of_links as f64)*BOLTZMANN_CONSTANT*temperature/6.0;
                    let residual_abs = &helmholtz_free_energy - &helmholtz_free_energy_from_initial_slope;
                    let residual_rel = &residual_abs/&helmholtz_free_energy;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let small_force = parameters.small*&force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let helmholtz_free_energy_per_link = model.helmholtz_free_energy_per_link(&small_force, temperature);
                    let helmholtz_free_energy_per_link_from_initial_slope = model.helmholtz_free_energy_per_link(&parameters.zero, temperature) + (small_force/BOLTZMANN_CONSTANT/temperature*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature/6.0;
                    let residual_abs = &helmholtz_free_energy_per_link - &helmholtz_free_energy_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&helmholtz_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let small_force = parameters.small*&force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_helmholtz_free_energy = model.relative_helmholtz_free_energy(&small_force, temperature);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let force = parameters.force_reference + parameters.force_scale*(0.5 - rng.gen::<f64>());
                    let small_force = parameters.small*&force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let relative_helmholtz_free_energy_per_link = model.relative_helmholtz_free_energy_per_link(&small_force, temperature);
                    let relative_helmholtz_free_energy_per_link_from_initial_slope = (small_force/BOLTZMANN_CONSTANT/temperature*link_length).powf(2.0)*BOLTZMANN_CONSTANT*temperature/6.0;
                    let residual_abs = &relative_helmholtz_free_energy_per_link - &relative_helmholtz_free_energy_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&relative_helmholtz_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_force = parameters.small*&nondimensional_force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_helmholtz_free_energy = model.nondimensional_helmholtz_free_energy(&small_nondimensional_force, temperature);
                    let nondimensional_helmholtz_free_energy_from_initial_slope = model.nondimensional_helmholtz_free_energy(&parameters.zero, temperature) + small_nondimensional_force.powf(2.0)*(number_of_links as f64)/6.0;
                    let residual_abs = &nondimensional_helmholtz_free_energy - &nondimensional_helmholtz_free_energy_from_initial_slope;
                    let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_force = parameters.small*&nondimensional_force;
                    let temperature = parameters.temperature_reference + parameters.temperature_scale*(0.5 - rng.gen::<f64>());
                    let nondimensional_helmholtz_free_energy_per_link = model.nondimensional_helmholtz_free_energy_per_link(&small_nondimensional_force, temperature);
                    let nondimensional_helmholtz_free_energy_per_link_from_initial_slope = model.nondimensional_helmholtz_free_energy_per_link(&parameters.zero, temperature) + small_nondimensional_force.powf(2.0)/6.0;
                    let residual_abs = &nondimensional_helmholtz_free_energy_per_link - &nondimensional_helmholtz_free_energy_per_link_from_initial_slope;
                    let residual_rel = &residual_abs/&nondimensional_helmholtz_free_energy_per_link;
                    assert!(residual_abs.abs() <= parameters.abs_tol_for_close);
                    assert!(residual_rel.abs() <= parameters.rel_tol_for_close);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_force = parameters.small*&nondimensional_force;
                    let nondimensional_relative_helmholtz_free_energy = model.nondimensional_relative_helmholtz_free_energy(&small_nondimensional_force);
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
                    let model = <$model>::init(number_of_links, link_length, hinge_mass);
                    let nondimensional_force = parameters.nondimensional_force_reference + parameters.nondimensional_force_scale*(0.5 - rng.gen::<f64>());
                    let small_nondimensional_force = parameters.small*&nondimensional_force;
                    let nondimensional_relative_helmholtz_free_energy_per_link = model.nondimensional_relative_helmholtz_free_energy_per_link(&small_nondimensional_force);
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
pub(crate) use isotensionalLegendre;