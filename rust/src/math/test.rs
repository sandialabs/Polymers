#![cfg(test)]

use rand::Rng;
use crate::math::integrate;

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
