#![cfg(test)]

use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};

pub struct Parameters
{
    pub number_of_loops: u16,
    pub abs_tol: f64,
    pub rel_tol: f64,
    pub force_scale: f64,
    pub nondimensional_force_scale: f64,
    pub temperature_scale: f64,
    pub hinge_mass_scale: f64,
    pub link_length_scale: f64,
    pub default_link_length: f64,
    pub minimum_number_of_links: u16,
    pub maximum_number_of_links: u16,
    pub default_number_of_links: u16,
    pub default_hinge_mass: f64,
}

impl Default for Parameters
{
    fn default() -> Self
    {
        Self
        {
            number_of_loops: 8,
            abs_tol: 1e-9,
            rel_tol: 1e-9,
            force_scale: 1e1,
            nondimensional_force_scale: 1e1,
            temperature_scale: 3e2,
            hinge_mass_scale: 1e0,
            link_length_scale: 1e0,
            default_link_length: 1e0,
            minimum_number_of_links: 8,
            maximum_number_of_links: 88,
            default_number_of_links: 8,
            default_hinge_mass: 1.0
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
                let link_length = parameters.link_length_scale*rng.gen::<f64>();
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
                let hinge_mass = parameters.hinge_mass_scale*rng.gen::<f64>();
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
                let hinge_mass = parameters.hinge_mass_scale*rng.gen::<f64>();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    for _ in 0..parameters.number_of_loops
                    {
                        let link_length = rng.gen::<f64>();
                        assert_eq!(link_length, <$model>::init(number_of_links, link_length, hinge_mass).link_length);
                    }
                }
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
                    for _ in 0..parameters.number_of_loops
                    {
                        let link_length = parameters.link_length_scale*rng.gen::<f64>();
                        let model = <$model>::init(number_of_links, link_length, parameters.default_hinge_mass);
                        assert_eq!(number_of_links + 1, model.random_configuration().len() as u16);
                    }
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
                    for _ in 0..parameters.number_of_loops
                    {
                        let link_length = parameters.link_length_scale*rng.gen::<f64>();
                        let force = parameters.force_scale*rng.gen::<f64>();
                        let temperature = parameters.temperature_scale*rng.gen::<f64>();
                        let model = <$model>::init(number_of_links, link_length, parameters.default_hinge_mass);
                        let end_to_end_length = model.end_to_end_length(&force, temperature);
                        let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                        let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_force);
                        let residual_abs = &end_to_end_length/link_length - &nondimensional_end_to_end_length;
                        let residual_rel = &residual_abs/&nondimensional_end_to_end_length;
                        assert!(residual_abs.abs() <= parameters.abs_tol);
                        assert!(residual_rel.abs() <= parameters.rel_tol);
                    }
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
                    for _ in 0..parameters.number_of_loops
                    {
                        let link_length = parameters.link_length_scale*rng.gen::<f64>();
                        let force = parameters.force_scale*rng.gen::<f64>();
                        let temperature = parameters.temperature_scale*rng.gen::<f64>();
                        let model = <$model>::init(number_of_links, link_length, parameters.default_hinge_mass);
                        let end_to_end_length_per_link = model.end_to_end_length_per_link(&force, temperature);
                        let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                        let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                        let residual_abs = &end_to_end_length_per_link/link_length - &nondimensional_end_to_end_length_per_link;
                        let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link;
                        assert!(residual_abs.abs() <= parameters.abs_tol);
                        assert!(residual_rel.abs() <= parameters.rel_tol);
                    }
                }
            }
            // #[test]
            // fn gibbs_free_energy()
            // {}
            #[test]
            fn relative_gibbs_free_energy()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    for _ in 0..parameters.number_of_loops
                    {
                        let link_length = parameters.link_length_scale*rng.gen::<f64>();
                        let force = parameters.force_scale*rng.gen::<f64>();
                        let temperature = parameters.temperature_scale*rng.gen::<f64>();
                        let model = <$model>::init(number_of_links, link_length, parameters.default_hinge_mass);
                        let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&force, temperature);
                        let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                        let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
                        let residual_abs = &relative_gibbs_free_energy/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_gibbs_free_energy;
                        let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy;
                        assert!(residual_abs.abs() <= parameters.abs_tol);
                        assert!(residual_rel.abs() <= parameters.rel_tol);
                    }
                }
            }
            // #[test]
            // fn gibbs_free_energy_per_link()
            // {}
            #[test]
            fn relative_gibbs_free_energy_per_link()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    for _ in 0..parameters.number_of_loops
                    {
                        let link_length = parameters.link_length_scale*rng.gen::<f64>();
                        let force = parameters.force_scale*rng.gen::<f64>();
                        let temperature = parameters.temperature_scale*rng.gen::<f64>();
                        let model = <$model>::init(number_of_links, link_length, parameters.default_hinge_mass);
                        let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&force, temperature);
                        let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                        let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
                        let residual_abs = &relative_gibbs_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_gibbs_free_energy_per_link;
                        let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
                        assert!(residual_abs.abs() <= parameters.abs_tol);
                        assert!(residual_rel.abs() <= parameters.rel_tol);
                    }
                }
            }
            mod legendre
            {
                use super::*;
                use rand::Rng;
                use crate::physics::single_chain::test::Parameters;
                // #[test]
                // fn helmholtz_free_energy()
                // {}
                #[test]
                fn relative_helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        for _ in 0..parameters.number_of_loops
                        {
                            let link_length = parameters.link_length_scale*rng.gen::<f64>();
                            let force = parameters.force_scale*rng.gen::<f64>();
                            let temperature = parameters.temperature_scale*rng.gen::<f64>();
                            let model = <$model>::init(number_of_links, link_length, parameters.default_hinge_mass);
                            let relative_helmholtz_free_energy = model.legendre.relative_helmholtz_free_energy(&force, temperature);
                            let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                            let nondimensional_relative_helmholtz_free_energy = model.legendre.nondimensional_relative_helmholtz_free_energy(&nondimensional_force);
                            let residual_abs = &relative_helmholtz_free_energy/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_helmholtz_free_energy;
                            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy;
                            assert!(residual_abs.abs() <= parameters.abs_tol);
                            assert!(residual_rel.abs() <= parameters.rel_tol);
                        }
                    }
                }
                // #[test]
                // fn helmholtz_free_energy_per_link()
                // {}
                #[test]
                fn relative_helmholtz_free_energy_per_link()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        for _ in 0..parameters.number_of_loops
                        {
                            let link_length = parameters.link_length_scale*rng.gen::<f64>();
                            let force = parameters.force_scale*rng.gen::<f64>();
                            let temperature = parameters.temperature_scale*rng.gen::<f64>();
                            let model = <$model>::init(number_of_links, link_length, parameters.default_hinge_mass);
                            let relative_helmholtz_free_energy_per_link = model.legendre.relative_helmholtz_free_energy_per_link(&force, temperature);
                            let nondimensional_force = force/BOLTZMANN_CONSTANT/temperature*link_length;
                            let nondimensional_relative_helmholtz_free_energy_per_link = model.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force);
                            let residual_abs = &relative_helmholtz_free_energy_per_link/BOLTZMANN_CONSTANT/temperature - &nondimensional_relative_helmholtz_free_energy_per_link;
                            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
                            assert!(residual_abs.abs() <= parameters.abs_tol);
                            assert!(residual_rel.abs() <= parameters.rel_tol);
                        }
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
                    for _ in 0..parameters.number_of_loops
                    {
                        let link_length = parameters.link_length_scale*rng.gen::<f64>();
                        let force = parameters.force_scale*rng.gen::<f64>();
                        let temperature = parameters.temperature_scale*rng.gen::<f64>();
                        let model = <$model>::init(number_of_links, link_length, parameters.default_hinge_mass);
                        let end_to_end_length = model.end_to_end_length(&force, temperature);
                        let end_to_end_length_per_link = model.end_to_end_length_per_link(&force, temperature);
                        let residual_abs = &end_to_end_length/(model.number_of_links as f64) - &end_to_end_length_per_link;
                        let residual_rel = &residual_abs/&end_to_end_length_per_link;
                        assert!(residual_abs.abs() <= parameters.abs_tol);
                        assert!(residual_rel.abs() <= parameters.rel_tol);
                    }
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
                    for _ in 0..parameters.number_of_loops
                    {
                        let link_length = parameters.link_length_scale*rng.gen::<f64>();
                        let nondimensional_force = parameters.nondimensional_force_scale*rng.gen::<f64>();
                        let model = <$model>::init(number_of_links, link_length, parameters.default_hinge_mass);
                        let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_force);
                        let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                        let residual_abs = &nondimensional_end_to_end_length/(model.number_of_links as f64) - &nondimensional_end_to_end_length_per_link;
                        let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link;
                        assert!(residual_abs.abs() <= parameters.abs_tol);
                        assert!(residual_rel.abs() <= parameters.rel_tol);
                    }
                }
            }
            // #[test]
            // fn gibbs_free_energy()
            // {}
            #[test]
            fn relative_gibbs_free_energy()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    for _ in 0..parameters.number_of_loops
                    {
                        let link_length = parameters.link_length_scale*rng.gen::<f64>();
                        let force = parameters.force_scale*rng.gen::<f64>();
                        let temperature = parameters.temperature_scale*rng.gen::<f64>();
                        let model = <$model>::init(number_of_links, link_length, parameters.default_hinge_mass);
                        let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&force, temperature);
                        let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&force, temperature);
                        let residual_abs = &relative_gibbs_free_energy/(model.number_of_links as f64) - &relative_gibbs_free_energy_per_link;
                        let residual_rel = &residual_abs/&relative_gibbs_free_energy_per_link;
                        assert!(residual_abs.abs() <= parameters.abs_tol);
                        assert!(residual_rel.abs() <= parameters.rel_tol);
                    }
                }
            }
            // #[test]
            // fn nondimensional_gibbs_free_energy()
            // {}
            #[test]
            fn nondimensional_relative_gibbs_free_energy()
            {
                let parameters = Parameters::default();
                let mut rng = rand::thread_rng();
                for _ in 0..parameters.number_of_loops
                {
                    let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                    for _ in 0..parameters.number_of_loops
                    {
                        let link_length = parameters.link_length_scale*rng.gen::<f64>();
                        let nondimensional_force = parameters.nondimensional_force_scale*rng.gen::<f64>();
                        let model = <$model>::init(number_of_links, link_length, parameters.default_hinge_mass);
                        let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
                        let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
                        let residual_abs = &nondimensional_relative_gibbs_free_energy/(model.number_of_links as f64) - &nondimensional_relative_gibbs_free_energy_per_link;
                        let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
                        assert!(residual_abs.abs() <= parameters.abs_tol);
                        assert!(residual_rel.abs() <= parameters.rel_tol);
                    }
                }
            }
            mod legendre
            {
                use super::*;
                use rand::Rng;
                use crate::physics::single_chain::test::Parameters;
                // #[test]
                // fn helmholtz_free_energy()
                // {}
                #[test]
                fn relative_helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        for _ in 0..parameters.number_of_loops
                        {
                            let link_length = parameters.link_length_scale*rng.gen::<f64>();
                            let force = parameters.force_scale*rng.gen::<f64>();
                            let temperature = parameters.temperature_scale*rng.gen::<f64>();
                            let model = <$model>::init(number_of_links, link_length, parameters.default_hinge_mass);
                            let relative_helmholtz_free_energy = model.legendre.relative_helmholtz_free_energy(&force, temperature);
                            let relative_helmholtz_free_energy_per_link = model.legendre.relative_helmholtz_free_energy_per_link(&force, temperature);
                            let residual_abs = &relative_helmholtz_free_energy/(model.number_of_links as f64) - &relative_helmholtz_free_energy_per_link;
                            let residual_rel = &residual_abs/&relative_helmholtz_free_energy_per_link;
                            assert!(residual_abs.abs() <= parameters.abs_tol);
                            assert!(residual_rel.abs() <= parameters.rel_tol);
                        }
                    }
                }
                // #[test]
                // fn nondimensional_helmholtz_free_energy()
                // {}
                #[test]
                fn nondimensional_relative_helmholtz_free_energy()
                {
                    let parameters = Parameters::default();
                    let mut rng = rand::thread_rng();
                    for _ in 0..parameters.number_of_loops
                    {
                        let number_of_links: u16 = rng.gen_range(parameters.minimum_number_of_links..parameters.maximum_number_of_links);
                        for _ in 0..parameters.number_of_loops
                        {
                            let link_length = parameters.link_length_scale*rng.gen::<f64>();
                            let nondimensional_force = parameters.nondimensional_force_scale*rng.gen::<f64>();
                            let model = <$model>::init(number_of_links, link_length, parameters.default_hinge_mass);
                            let nondimensional_relative_helmholtz_free_energy = model.legendre.nondimensional_relative_helmholtz_free_energy(&nondimensional_force);
                            let nondimensional_relative_helmholtz_free_energy_per_link = model.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force);
                            let residual_abs = &nondimensional_relative_helmholtz_free_energy/(model.number_of_links as f64) - &nondimensional_relative_helmholtz_free_energy_per_link;
                            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
                            assert!(residual_abs.abs() <= parameters.abs_tol);
                            assert!(residual_rel.abs() <= parameters.rel_tol);
                        }
                    }
                }
            }
        }
        mod relative
        {
            // add in tests at zero force for absolute free energy / partition functions (zero force == free ensemble, dont need to make that one)
            // going to need hinge mass and Planck constant in order to nondimensionalize things
            // at that point it might be good to use Boltzmann as well (put in physics) and use temperature rather than inverse temperature...
            mod legendre
            {}
        }
    }
}
pub(crate) use isotensional;
