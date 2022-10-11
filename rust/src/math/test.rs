#![cfg(test)]

use crate::math::factorial;

#[test]
fn largest_factorial()
{
    factorial(34);
}

mod inverse_langevin
{
    mod scalar
    {
        use rand::Rng;
        use crate::math::langevin;
        use crate::math::inverse_langevin;

        #[test]
        fn inverse()
        {
            let mut rng = rand::thread_rng();
            for _ in 0..88888
            {
                let y = rng.gen::<f64>();
                let x = inverse_langevin(&y, 1e-4);
                let f = langevin(&x);
                let residual_abs = &y - &f;
                let residual_rel = &residual_abs/&y;
                assert!(residual_abs.abs() <= 1e-4);
                assert!(residual_rel.abs() <= 1e-4);
            }
        }
    }
}

mod approximate_inverse_langevin
{
    mod scalar
    {
        use rand::Rng;
        use crate::math::langevin;
        use crate::math::approximate_inverse_langevin;

        #[test]
        fn inverse()
        {
            let mut rng = rand::thread_rng();
            for _ in 0..88888
            {
                let y = rng.gen::<f64>();
                let x = approximate_inverse_langevin(&y);
                let f = langevin(&x);
                let residual_abs = &y - &f;
                let residual_rel = &residual_abs/&y;
                assert!(residual_abs.abs() <= 1e-3);
                assert!(residual_rel.abs() <= 1e-3);
            }
        }
    }
    mod array_1d
    {
        use ndarray::Array;
        use crate::math::langevin;
        use crate::math::approximate_inverse_langevin;

        #[test]
        fn inverse()
        {
            let y = Array::<f64, _>::linspace(1e-6, 1.0, 88888);
            let x = approximate_inverse_langevin(&y);
            let f = langevin(&x);
            let residual_abs = &y - &f;
            let residual_rel = &residual_abs/&y;
            for residual_abs_i in residual_abs.iter()
            {
                assert!(residual_abs_i.abs() <= 1e-3);
            }
            for residual_rel_i in residual_rel.iter()
            {
                assert!(residual_rel_i.abs() <= 1e-3);
            }
        }
    }
}

mod langevin
{

    mod scalar
    {
        use rand::Rng;
        use crate::math::langevin;

        static SIZE: u32 = 88888;
        static SCALE: f64 = 1e-3;
        static ABS_TOL: f64 = 1e-7;
        static REL_TOL: f64 = 1e-5;

        #[test]
        fn scalar_small()
        {
            for _ in 1..SIZE
            {
                let random_float_0_1: f64 = rand::thread_rng().gen();
                let x = SCALE*(1.0 + 9.0*random_float_0_1);
                let y = langevin(&x);
                let residual_abs = &y - &x/3.0;
                let residual_rel = &residual_abs/&y;
                assert!(residual_abs.abs() <= ABS_TOL);
                assert!(residual_rel.abs() <= REL_TOL);
            }
        }
    }

    mod array_1d
    {
        use ndarray::Array;
        use crate::math::langevin;

        static SIZE: usize = 88888;
        static SCALE: f64 = 1e-3;
        static ABS_TOL: f64 = 1e-7;
        static REL_TOL: f64 = 1e-5;

        #[test]
        fn array_small()
        {
            let x = SCALE*Array::<f64, _>::linspace(1.0, 10.0, SIZE);
            let y = langevin(&x);
            let residual_abs = &y - &x/3.0;
            let residual_rel = &residual_abs/&y;
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

    mod array_2d
    {
        use ndarray::Array;
        use crate::math::langevin;

        static SIZE: usize = 88888;
        static QUARTER_SIZE: usize = 22222;
        static SCALE: f64 = 1e-3;
        static ABS_TOL: f64 = 1e-7;
        static REL_TOL: f64 = 1e-5;

        #[test]
        fn array_small()
        {
            let x = SCALE*Array::<f64, _>::linspace(1.0, 10.0, SIZE).into_shape((QUARTER_SIZE, 4)).unwrap();
            let y = langevin(&x);
            let residual_abs = &y - &x/3.0;
            let residual_rel = &residual_abs/&y;
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
