#![cfg(test)]

use rand::Rng;
use crate::math::integrate;
use crate::math::invert;
use crate::math::sequential_uniform_random;

#[test]
fn test_sequential_uniform_random()
{
    let mut x: f64 = 1.0;
    let mut x_bar = 0.0;
    let mut x_squared_bar = 0.0;
    let num = 888888;
    for _ in 0..num
    {
        x = sequential_uniform_random(x);
        x_bar += x/(num as f64);
        x_squared_bar += x*x/(num as f64);
        assert!(x >= 0.0 && x <= 1.0);
    }
    assert!((x_bar - 0.5).abs() <= 1e-3);
    assert!((x_squared_bar - 1.0/3.0).abs() <= 1e-3);
}

#[test]
fn integrate_cosh()
{
    fn cosh(x: f64) -> f64 {x.cosh()}
    let mut rng = rand::thread_rng();
    for _ in 0..88
    {
        let limit = rng.gen::<f64>();
        let integral = integrate(cosh, &-limit, &limit, &10000);
        assert!((integral - 2.0*limit.sinh()).abs() <= 1e-8);
    }
}

#[test]
fn invert_sinh()
{
    fn sinh(x: f64) -> f64 {x.sinh()}
    let mut rng = rand::thread_rng();
    for _ in 0..8888
    {
        let y = rng.gen::<f64>();
        let x = invert(y, sinh, 0.0);
        let residual_abs = &y - &x.sinh();
        let residual_rel = &residual_abs/&y;
        assert!(residual_abs.abs() <= 1e-4);
        assert!(residual_rel.abs() <= 1e-4);
    }
}

mod inverse_langevin
{
    use rand::Rng;
    use crate::math::inverse_langevin;
    #[test]
    #[ignore]
    fn inverse()
    {
        let mut rng = rand::thread_rng();
        for _ in 0..8888
        {
            let y = rng.gen::<f64>();
            let x = inverse_langevin(&y, 1e-4);
            let f = 1.0/x.tanh() - 1.0/x;
            let residual_abs = &y - &f;
            let residual_rel = &residual_abs/&y;
            assert!(residual_abs.abs() <= 1e-4);
            assert!(residual_rel.abs() <= 1e-4);
        }
    }
}

mod approximate_inverse_langevin
{
    use rand::Rng;
    use crate::math::approximate_inverse_langevin;
    #[test]
    fn inverse()
    {
        let mut rng = rand::thread_rng();
        for _ in 0..8888
        {
            let y = rng.gen::<f64>();
            let x = approximate_inverse_langevin(&y);
            let f = 1.0/x.tanh() - 1.0/x;
            let residual_abs = &y - &f;
            let residual_rel = &residual_abs/&y;
            assert!(residual_abs.abs() <= 1e-3);
            assert!(residual_rel.abs() <= 1e-3);
        }
    }
}
