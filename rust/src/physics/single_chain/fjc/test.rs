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
mod implementations
{
    use super::*;
    mod thermodynamics
    {
        use super::*;
        mod isometric
        {
            use super::*;
            #[test]
            fn access()
            {
                let parameters = Parameters::default();
                let _ = FJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference).thermodynamics.isometric.nondimensional_force(&parameters.nondimensional_end_to_end_length_per_link_reference);
            }
            mod legendre
            {
                use super::*;
                #[test]
                fn access()
                {
                    let parameters = Parameters::default();
                    let _ = FJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference).thermodynamics.isometric.legendre.nondimensional_relative_gibbs_free_energy_per_link(&parameters.nondimensional_end_to_end_length_per_link_reference);
                }
            }
        }
        mod isotensional
        {
            use super::*;
            #[test]
            fn access()
            {
                let parameters = Parameters::default();
                let _ = FJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference).thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(&parameters.nondimensional_force_reference);
            }
            mod legendre
            {
                use super::*;
                #[test]
                fn access()
                {
                    let parameters = Parameters::default();
                    let _ = FJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference).thermodynamics.isotensional.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&parameters.nondimensional_force_reference);
                }
            }
        }
        mod modified_canonical
        {
            use super::*;
            #[test]
            fn access()
            {
                let parameters = Parameters::default();
                let _ = FJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference).thermodynamics.modified_canonical.nondimensional_force(&parameters.nondimensional_potential_distance_small, &parameters.nondimensional_potential_stiffness_reference);
            }
            mod asymptotic
            {
                use super::*;
                mod strong_potential
                {
                    use super::*;
                    #[test]
                    fn access()
                    {
                        let parameters = Parameters::default();
                        let _ = FJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference).thermodynamics.modified_canonical.asymptotic.strong_potential.nondimensional_force(&parameters.nondimensional_potential_distance_small, &parameters.nondimensional_potential_stiffness_reference);
                    }
                }
                mod weak_potential
                {
                    use super::*;
                    #[test]
                    fn access()
                    {
                        let parameters = Parameters::default();
                        let _ = FJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference).thermodynamics.modified_canonical.asymptotic.weak_potential.nondimensional_end_to_end_length_per_link(&parameters.nondimensional_potential_distance_reference, &parameters.nondimensional_potential_stiffness_reference);
                    }
                }
            }
        }
    }
}