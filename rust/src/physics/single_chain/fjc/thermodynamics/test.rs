#![cfg(test)]

mod init {

    use rand::prelude::*;
    use crate::physics::single_chain::fjc::thermodynamics::Thermodynamics;

    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        for _ in 0..88
        {
            let number_of_links: u16 = rng.gen_range(8..88);
            assert_eq!(number_of_links, Thermodynamics::init(number_of_links).number_of_links);
        }
    }

    #[test]
    fn isometric()
    {
        let _ = Thermodynamics::init(8).isometric;
    }

    #[test]
    fn isotensional()
    {
        let _ = Thermodynamics::init(8).isotensional;
    }
}

mod legendre_transformation {

    use rand::prelude::*;
    use crate::physics::single_chain::fjc::thermodynamics::Thermodynamics;

    #[test]
    fn helmholtz_free_energy()
    {
        let mut rng = rand::thread_rng();
        for _ in 0..88
        {
            let number_of_links: u16 = rng.gen_range(8..88);
            let nondimensional_force: f64 = rng.gen();
            let model = Thermodynamics::init(number_of_links);
            assert_eq!(
                model.isometric.nondimensional_relative_helmholtz_free_energy_per_link_legendre_transformation(nondimensional_force),
                model.isotensional.nondimensional_relative_gibbs_free_energy_per_link(nondimensional_force)
                    + nondimensional_force*model.isotensional.nondimensional_end_to_end_length_per_link(nondimensional_force)
            );
        }
    }
}