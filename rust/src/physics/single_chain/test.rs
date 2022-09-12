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
                assert_eq!(number_of_links, <$model>::init(number_of_links).number_of_links);
            }
        }
    }	
}
pub(crate) use base;

macro_rules! single_chain
{
    ( $model:ty ) =>
    {
        #[test]
        fn mechanics()
        {
            let _ = <$model>::init(8).mechanics;
        }
        #[test]
        fn thermodynamics()
        {
            let _ = <$model>::init(8).thermodynamics;
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
                let model = <$model>::init(number_of_links);
                assert_eq!(number_of_links + 1, model.random_configuration().len() as u16);
            }
        }
    }
}
pub(crate) use mechanics;

macro_rules! thermodynamics
{
    ( $model:ty ) =>
    {
        #[test]
        fn isometric()
        {
            let _ = <$model>::init(8).isometric;
        }
        #[test]
        fn isotensional()
        {
            let _ = <$model>::init(8).isotensional;
        }
        #[test]
        fn legendre_transformation_nondimensional_relative_helmholtz_free_energy_per_link()
        {
            let mut rng = rand::thread_rng();
            for _ in 0..8
            {
                let number_of_links: u16 = rng.gen_range(8..88);
                let model = <$model>::init(number_of_links);
                model.isotensional.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&1.0);
                // test whether approximation becomes accurate for large number_of_links
            }
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
    ( $model:ty ) => {}
}
pub(crate) use isotensional;
