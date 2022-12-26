#![cfg(test)]
use super::*;
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
            number_of_loops: 888,
            abs_tol: 1e-8,
            rel_tol: 1e-6,
            log_log_tol: 5e-2,
            log_log_scale: 12e-1,
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
            nondimensional_force_reference: 5e1,
            nondimensional_force_scale: 1e2,
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
        let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference);
    }
    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            assert_eq!(number_of_links, EFJC::init(number_of_links, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).number_of_links);
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
            assert_eq!(link_length, EFJC::init(parameters.number_of_links_minimum, link_length, parameters.hinge_mass_reference, parameters.link_stiffness_reference).link_length);
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
            assert_eq!(hinge_mass, EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, hinge_mass, parameters.link_stiffness_reference).hinge_mass);
        }
    }
    #[test]
    fn link_stiffness()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(link_stiffness, EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, link_stiffness).link_stiffness);
        }
    }
    #[test]
    fn number_of_links_and_link_length_and_hinge_mass_and_link_stiffness()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_length = rng.gen::<f64>();
            let link_stiffness = parameters.link_stiffness_reference + parameters.link_stiffness_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(link_length, EFJC::init(number_of_links, link_length, hinge_mass, link_stiffness).link_length);
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
            #[test]
            fn access()
            {
                let parameters = Parameters::default();
                let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).thermodynamics.isotensional.nondimensional_end_to_end_length_per_link(&parameters.nondimensional_force_reference, &parameters.temperature_reference);
            }
            mod legendre
            {
                use super::*;
                #[test]
                fn access()
                {
                    let parameters = Parameters::default();
                    let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).thermodynamics.isotensional.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&parameters.nondimensional_force_reference, &parameters.temperature_reference);
                }
            }
            use super::*;
            mod asymptotic
            {
                use super::*;
                #[test]
                fn access()
                {
                    let parameters = Parameters::default();
                    let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).thermodynamics.isotensional.asymptotic.nondimensional_end_to_end_length_per_link(&parameters.nondimensional_force_reference, &parameters.temperature_reference);
                }
                mod legendre
                {
                    use super::*;
                    #[test]
                    fn access()
                    {
                        let parameters = Parameters::default();
                        let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).thermodynamics.isotensional.asymptotic.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&parameters.nondimensional_force_reference, &parameters.temperature_reference);
                    }
                }
                mod alternative
                {
                    use super::*;
                    #[test]
                    fn access()
                    {
                        let parameters = Parameters::default();
                        let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).thermodynamics.isotensional.asymptotic.alternative.nondimensional_end_to_end_length_per_link(&parameters.nondimensional_force_reference, &parameters.temperature_reference);
                    }
                    mod legendre
                    {
                        use super::*;
                        #[test]
                        fn access()
                        {
                            let parameters = Parameters::default();
                            let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).thermodynamics.isotensional.asymptotic.alternative.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&parameters.nondimensional_force_reference, &parameters.temperature_reference);
                        }
                    }
                }
                mod reduced
                {
                    use super::*;
                    #[test]
                    fn access()
                    {
                        let parameters = Parameters::default();
                        let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).thermodynamics.isotensional.asymptotic.reduced.nondimensional_end_to_end_length_per_link(&parameters.nondimensional_force_reference, &parameters.temperature_reference);
                    }
                    mod legendre
                    {
                        use super::*;
                        #[test]
                        fn access()
                        {
                            let parameters = Parameters::default();
                            let _ = EFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.link_stiffness_reference).thermodynamics.isotensional.asymptotic.reduced.legendre.nondimensional_relative_helmholtz_free_energy_per_link(&parameters.nondimensional_force_reference, &parameters.temperature_reference);
                        }
                    }
                }
            }
        }
    }
}
