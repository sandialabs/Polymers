#![cfg(test)]

mod init {

    use rand::prelude::*;
    use crate::physics::single_chain::fjc::thermodynamics::isotensional::Isotensional;

    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        for _ in 0..8
        {
            let number_of_links: u16 = rng.gen_range(8..88);
            assert_eq!(number_of_links, Isotensional::init(number_of_links).number_of_links);
        }
    }

}
mod helmholtz_free_energy {

    mod per_link {

        use ndarray::Array;
        use rand::prelude::*;
        use crate::physics::single_chain::fjc::thermodynamics::isotensional::Isotensional;

        static SCALE: f64 = 1e1;
        static ABS_TOL: f64 = 1e-7;
        static REL_TOL: f64 = 1e-5;

        #[test]
        fn scalar()
        {
            let mut rng = rand::thread_rng();
            for _ in 0..8
            {
                let number_of_links: u16 = rng.gen_range(8..88);
                let nondimensional_force = SCALE*rng.gen::<f64>();
                let model = Isotensional::init(number_of_links);
                let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
                let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
                let residual_abs = &nondimensional_relative_gibbs_free_energy/(model.number_of_links as f64) - &nondimensional_relative_gibbs_free_energy_per_link;
                let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
                assert!(residual_abs.abs() <= ABS_TOL);
                assert!(residual_rel.abs() <= REL_TOL);
            }
        }

        #[test]
        fn array_1d()
        {
            let mut rng = rand::thread_rng();
            for _ in 0..8
            {
                let number_of_links: u16 = rng.gen_range(8..88);
                let nondimensional_force = SCALE*Array::<f64, _>::linspace(1.0, 10.0, 88);
                let model = Isotensional::init(number_of_links);
                let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
                let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
                let residual_abs = &nondimensional_relative_gibbs_free_energy/(model.number_of_links as f64) - &nondimensional_relative_gibbs_free_energy_per_link;
                let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
                for residual_abs_i in residual_abs.iter()
                {
                    assert!(residual_abs_i.abs() <= ABS_TOL);
                }
                for residual_rel_i in residual_rel.iter()
                {
                    assert!(residual_rel_i.abs() <= REL_TOL);
                }
            }
        }

        #[test]
        fn array_2d()
        {
            let mut rng = rand::thread_rng();
            for _ in 0..8
            {
                let number_of_links: u16 = rng.gen_range(8..88);
                let nondimensional_force = SCALE*Array::<f64, _>::linspace(1.0, 10.0, 88).into_shape((22, 4)).unwrap();
                let model = Isotensional::init(number_of_links);
                let nondimensional_relative_gibbs_free_energy = model.nondimensional_relative_gibbs_free_energy(&nondimensional_force);
                let nondimensional_relative_gibbs_free_energy_per_link = model.nondimensional_relative_gibbs_free_energy_per_link(&nondimensional_force);
                let residual_abs = &nondimensional_relative_gibbs_free_energy/(model.number_of_links as f64) - &nondimensional_relative_gibbs_free_energy_per_link;
                let residual_rel = &residual_abs/&nondimensional_relative_gibbs_free_energy_per_link;
                for residual_abs_i in residual_abs.iter()
                {
                    assert!(residual_abs_i.abs() <= ABS_TOL);
                }
                for residual_rel_i in residual_rel.iter()
                {
                    assert!(residual_rel_i.abs() <= REL_TOL);
                }
            }
        }
    }
}
