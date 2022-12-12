#![cfg(test)]
use super::*;
use crate::physics::single_chain::test::Parameters;
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
mod implementations
{
    use super::*;
    mod thermodynamics
    {
        use super::*;
        mod isotensional
        {
            use super::*;
            #[test]
            fn access()
            {
                let parameters = Parameters::default();
                let _ = SWFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.well_width_reference).thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(&parameters.nondimensional_force_reference);
            }
            mod legendre
            {
                use super::*;
                #[test]
                fn access()
                {
                    let parameters = Parameters::default();
                    let _ = SWFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.well_width_reference).thermodynamics.isotensional.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&parameters.nondimensional_force_reference);
                }
            }
        }
    }
}