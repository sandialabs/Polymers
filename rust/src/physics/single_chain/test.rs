#![cfg(test)]

mod langevin {

    mod scalar {

        use rand::prelude::*;
        use crate::physics::single_chain::Math;

        static SIZE: u32 = 88888;
        static SCALE: f64 = 1e-3;
        static ABS_TOL: f64 = 1e-7;
        static REL_TOL: f64 = 1e-5;

        #[test]
        fn scalar_small() {
            for _ in 1..SIZE {
                let random_float_0_1: f64 = rand::thread_rng().gen();
                let eta = SCALE*(1.0 + 9.0*random_float_0_1);
                let y = Math::langevin(&eta);
                let residual_abs = &y - &eta/3.0;
                let residual_rel = &residual_abs/&y;
                assert!(residual_abs.abs() <= ABS_TOL);
                assert!(residual_rel.abs() <= REL_TOL);
            }
        }

    }

    mod array {

        use ndarray::Array;
        use crate::physics::single_chain::Math;

        static SIZE: usize = 88888;
        static SCALE: f64 = 1e-3;
        static ABS_TOL: f64 = 1e-7;
        static REL_TOL: f64 = 1e-5;

        #[test]
        fn array_small() {
            let eta = SCALE*Array::<f64, _>::linspace(1.0, 10.0, SIZE).into_shape((SIZE, 1)).unwrap();
            let y = Math::langevin(&eta);
            let residual_abs = &y - &eta/3.0;
            let residual_rel = &residual_abs/&y;
            for residual_abs_i in residual_abs.iter() {
                assert!(residual_abs_i.abs() <= ABS_TOL);
            }
            for residual_rel_i in residual_rel.iter() {
                assert!(residual_rel_i.abs() <= REL_TOL);
            }
        }

    }

}
