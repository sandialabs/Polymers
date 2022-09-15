macro_rules! base
{
    ( $model:ty ) =>
    {
        use rand::prelude::*;
        #[test]
        fn number_of_links()
        {
            let mut rng = rand::thread_rng();
            for _ in 0..8
            {
                let number_of_links: u16 = rng.gen_range(8..88);
                assert_eq!(number_of_links, <$model>::init(number_of_links, 1.0).number_of_links);
            }
        }
        #[test]
        fn link_length()
        {
            let mut rng = rand::thread_rng();
            for _ in 0..8
            {
                let link_length = rng.gen::<f64>();
                assert_eq!(link_length, <$model>::init(8, link_length).link_length);
            }
        }
        #[test]
        fn number_of_links_and_link_length()
        {
            let mut rng = rand::thread_rng();
            for _ in 0..8
            {
                let number_of_links: u16 = rng.gen_range(8..88);
                for _ in 0..8
                {
                    let link_length = rng.gen::<f64>();
                    assert_eq!(link_length, <$model>::init(number_of_links, link_length).link_length);
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
        mod init{
            use super::*;
            #[test]
            fn mechanics()
            {
                let _ = <$model>::init(8, 1.0).mechanics;
            }
            #[test]
            fn thermodynamics()
            {
                let _ = <$model>::init(8, 1.0).thermodynamics;
            }
        }
    }	
}
pub(crate) use single_chain;

macro_rules! mechanics
{
    ( $model:ty ) =>
    {
        #[test]
        fn number_of_configurations()
        {
            let mut rng = rand::thread_rng();
            for _ in 0..8
            {
                let number_of_links: u16 = rng.gen_range(8..88);
                for _ in 0..8
                {
                    let link_length = rng.gen::<f64>();
                    let model = <$model>::init(number_of_links, link_length);
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
        mod init{
            use super::*;
            #[test]
            fn isometric()
            {
                let _ = <$model>::init(8, 1.0).isometric;
            }
            #[test]
            fn isotensional()
            {
                let _ = <$model>::init(8, 1.0).isotensional;
            }
        }
        mod legendre_transformation{
            #[test]
            fn nondimensional_end_to_end_length()
            {}
            #[test]
            fn nondimensional_end_to_end_length_per_link()
            {}
            #[test]
            fn nondimensional_force()
            {}
            #[test]
            fn nondimensional_force_per_link()
            {}
            #[test]
            fn nondimensional_relative_helmholtz_free_energy()
            {}
            #[test]
            fn nondimensional_relative_helmholtz_free_energy_per_link()
            {}
            #[test]
            fn nondimensional_relative_gibbs_free_energy()
            {}
            #[test]
            fn nondimensional_relative_gibbs_free_energy_per_link()
            {}
        }
    }
}
pub(crate) use thermodynamics;

macro_rules! isometric
{
    ( $model:ty ) => {}
}
pub(crate) use isometric;

macro_rules! isotensional
{
    ( $model:ty ) =>
    {
        mod nondimensional
        {
            use super::*;
            use rand::prelude::*;
            static SCALE: f64 = 1e1;
            static ABS_TOL: f64 = 1e-7;
            static REL_TOL: f64 = 1e-5;
            #[test]
            fn end_to_end_length()
            {
                // let mut rng = rand::thread_rng();
                // for _ in 0..8
                // {
                //     let number_of_links: u16 = rng.gen_range(8..88);
                //     for _ in 0..8
                //     {
                //         let link_length = rng.gen::<f64>();
                //         let force = SCALE*rng.gen::<f64>();
                //         let inverse_temperature = SCALE*rng.gen::<f64>();
                //         let model = <$model>::init(number_of_links, link_length);
                //         let end_to_end_length = model.end_to_end_length(&force, inverse_temperature);
                //         let nondimensional_force = force*inverse_temperature*link_length;
                //         let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_force);
                //         let residual_abs = &end_to_end_length/link_length - &nondimensional_end_to_end_length;
                //         let residual_rel = &residual_abs/&nondimensional_end_to_end_length;
                //         assert!(residual_abs.abs() <= ABS_TOL);
                //         assert!(residual_rel.abs() <= REL_TOL);
                //     }
                // }
            }
            #[test]
            fn end_to_end_length_per_link()
            {
                let mut rng = rand::thread_rng();
                for _ in 0..8
                {
                    let number_of_links: u16 = rng.gen_range(8..88);
                    for _ in 0..8
                    {
                        let link_length = rng.gen::<f64>();
                        let force = SCALE*rng.gen::<f64>();
                        let inverse_temperature = SCALE*rng.gen::<f64>();
                        let model = <$model>::init(number_of_links, link_length);
                        let end_to_end_length_per_link = model.end_to_end_length_per_link(&force, inverse_temperature);
                        let nondimensional_force = force*inverse_temperature*link_length;
                        let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                        let residual_abs = &end_to_end_length_per_link/link_length - &nondimensional_end_to_end_length_per_link;
                        let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link;
                        assert!(residual_abs.abs() <= ABS_TOL);
                        assert!(residual_rel.abs() <= REL_TOL);
                    }
                }
            }
            // #[test]
            // fn gibbs_free_energy()
            // {

            // }
            // #[test]
            // fn relative_gibbs_free_energy()
            // {

            // }
            // #[test]
            // fn gibbs_free_energy_per_link()
            // {

            // }
            // #[test]
            // fn relative_gibbs_free_energy_per_link()
            // {

            // }
            mod legendre
            {
                // #[test]
                // fn helmholtz_free_energy()
                // {
    
                // }
                // #[test]
                // fn relative_helmholtz_free_energy()
                // {
    
                // }
                // #[test]
                // fn helmholtz_free_energy_per_link()
                // {
    
                // }
                // #[test]
                // fn relative_helmholtz_free_energy_per_link()
                // {
    
                // }
            }
        }
        mod per_link
        {
            use super::*;
            use rand::prelude::*;
            static SCALE: f64 = 1e1;
            static ABS_TOL: f64 = 1e-7;
            static REL_TOL: f64 = 1e-5;
            #[test]
            fn end_to_end_length()
            {
                let mut rng = rand::thread_rng();
                for _ in 0..8
                {
                    let number_of_links: u16 = rng.gen_range(8..88);
                    for _ in 0..8
                    {
                        let link_length = rng.gen::<f64>();
                        let force = SCALE*rng.gen::<f64>();
                        let inverse_temperature = SCALE*rng.gen::<f64>();
                        let model = <$model>::init(number_of_links, link_length);
                        let end_to_end_length = model.end_to_end_length(&force, inverse_temperature);
                        let end_to_end_length_per_link = model.end_to_end_length_per_link(&force, inverse_temperature);
                        let residual_abs = &end_to_end_length/(model.number_of_links as f64) - &end_to_end_length_per_link;
                        let residual_rel = &residual_abs/&end_to_end_length_per_link;
                        assert!(residual_abs.abs() <= ABS_TOL);
                        assert!(residual_rel.abs() <= REL_TOL);
                    }
                }
            }
            #[test]
            fn nondimensional_end_to_end_length()
            {
                let mut rng = rand::thread_rng();
                for _ in 0..8
                {
                    let number_of_links: u16 = rng.gen_range(8..88);
                    for _ in 0..8
                    {
                        let link_length = rng.gen::<f64>();
                        let nondimensional_force = SCALE*rng.gen::<f64>();
                        let model = <$model>::init(number_of_links, link_length);
                        let nondimensional_end_to_end_length = model.nondimensional_end_to_end_length(&nondimensional_force);
                        let nondimensional_end_to_end_length_per_link = model.nondimensional_end_to_end_length_per_link(&nondimensional_force);
                        let residual_abs = &nondimensional_end_to_end_length/(model.number_of_links as f64) - &nondimensional_end_to_end_length_per_link;
                        let residual_rel = &residual_abs/&nondimensional_end_to_end_length_per_link;
                        assert!(residual_abs.abs() <= ABS_TOL);
                        assert!(residual_rel.abs() <= REL_TOL);
                    }
                }
            }
            // #[test]
            // fn gibbs_free_energy()
            // {

            // }
            #[test]
            fn relative_gibbs_free_energy()
            {
                let mut rng = rand::thread_rng();
                for _ in 0..8
                {
                    let number_of_links: u16 = rng.gen_range(8..88);
                    for _ in 0..8
                    {
                        let link_length = rng.gen::<f64>();
                        let force = SCALE*rng.gen::<f64>();
                        let inverse_temperature = SCALE*rng.gen::<f64>();
                        let model = <$model>::init(number_of_links, link_length);
                        let relative_gibbs_free_energy = model.relative_gibbs_free_energy(&force, inverse_temperature);
                        let relative_gibbs_free_energy_per_link = model.relative_gibbs_free_energy_per_link(&force, inverse_temperature);
                        let residual_abs = &relative_gibbs_free_energy/(model.number_of_links as f64) - &relative_gibbs_free_energy_per_link;
                        let residual_rel = &residual_abs/&relative_gibbs_free_energy_per_link;
                        assert!(residual_abs.abs() <= ABS_TOL);
                        assert!(residual_rel.abs() <= REL_TOL);
                    }
                }
            }
            // #[test]
            // fn nondimensional_gibbs_free_energy()
            // {

            // }
            #[test]
            fn nondimensional_relative_gibbs_free_energy()
            {
                let mut rng = rand::thread_rng();
                for _ in 0..8
                {
                    let number_of_links: u16 = rng.gen_range(8..88);
                    for _ in 0..8
                    {
                        let link_length = rng.gen::<f64>();
                        let nondimensional_force = SCALE*rng.gen::<f64>();
                        let model = <$model>::init(number_of_links, link_length);
                        let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
                        let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
                        let residual_abs = &nondimensional_relative_gibbs_free_energy/(model.number_of_links as f64) - &nondimensional_relative_gibbs_free_energy_per_link;
                        let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
                        assert!(residual_abs.abs() <= ABS_TOL);
                        assert!(residual_rel.abs() <= REL_TOL);
                    }
                }
            }
            mod legendre
            {
                use super::*;
                use rand::prelude::*;
                static SCALE: f64 = 1e1;
                static ABS_TOL: f64 = 1e-7;
                static REL_TOL: f64 = 1e-5;
                // #[test]
                // fn helmholtz_free_energy()
                // {
    
                // }
                #[test]
                fn relative_helmholtz_free_energy()
                {
                    let mut rng = rand::thread_rng();
                    for _ in 0..8
                    {
                        let number_of_links: u16 = rng.gen_range(8..88);
                        for _ in 0..8
                        {
                            let link_length = rng.gen::<f64>();
                            let force = SCALE*rng.gen::<f64>();
                            let inverse_temperature = SCALE*rng.gen::<f64>();
                            let model = <$model>::init(number_of_links, link_length);
                            let relative_helmholtz_free_energy = model.legendre.relative_helmholtz_free_energy(&force, inverse_temperature);
                            let relative_helmholtz_free_energy_per_link = model.legendre.relative_helmholtz_free_energy_per_link(&force, inverse_temperature);
                            let residual_abs = &relative_helmholtz_free_energy/(model.number_of_links as f64) - &relative_helmholtz_free_energy_per_link;
                            let residual_rel = &residual_abs/&relative_helmholtz_free_energy_per_link;
                            assert!(residual_abs.abs() <= ABS_TOL);
                            assert!(residual_rel.abs() <= REL_TOL);
                        }
                    }
                }
                // #[test]
                // fn nondimensional_helmholtz_free_energy()
                // {
    
                // }
                #[test]
                fn nondimensional_relative_helmholtz_free_energy()
                {
                    let mut rng = rand::thread_rng();
                    for _ in 0..8
                    {
                        let number_of_links: u16 = rng.gen_range(8..88);
                        for _ in 0..8
                        {
                            let link_length = rng.gen::<f64>();
                            let nondimensional_force = SCALE*rng.gen::<f64>();
                            let model = <$model>::init(number_of_links, link_length);
                            let nondimensional_relative_helmholtz_free_energy = model.legendre.nondimensional_relative_helmholtz_free_energy(&nondimensional_force);
                            let nondimensional_relative_helmholtz_free_energy_per_link = model.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&nondimensional_force);
                            let residual_abs = &nondimensional_relative_helmholtz_free_energy/(model.number_of_links as f64) - &nondimensional_relative_helmholtz_free_energy_per_link;
                            let residual_rel = &residual_abs/&nondimensional_relative_helmholtz_free_energy_per_link;
                            assert!(residual_abs.abs() <= ABS_TOL);
                            assert!(residual_rel.abs() <= REL_TOL);
                        }
                    }
                }
            }
        }
        mod relative
        {
            // add in tests at zero force for absolute free energy / partition functions (zero force == free ensemble, dont need to make that one)
            mod legendre
            {

            }
        }
    }
}
pub(crate) use isotensional;