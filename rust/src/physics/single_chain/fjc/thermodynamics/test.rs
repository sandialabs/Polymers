#![cfg(test)]

mod init {

    use rand::prelude::*;
    use crate::physics::single_chain::fjc::thermodynamics::Thermodynamics;

    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        for _ in 0..8
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

mod legendre {

    use crate::physics::single_chain::fjc::thermodynamics::Thermodynamics;

    #[test]
    fn nondimensional_relative_helmholtz_free_energy_per_link()
    {
        let model = Thermodynamics::init(8);
        model.isotensional.legendre.nondimensional_relative_helmholtz_free_energy_per_link(1.0);
        // test whether approximation becomes accurate for large number_of_links
    }
}