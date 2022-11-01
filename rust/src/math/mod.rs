pub mod test;

pub fn invert<F>(value: f64, function: F, mut guess: f64) -> f64
where F: Fn(f64) -> f64
{
    let mut residual: f64 = 1.0;
    while residual.abs()/value > 1e-4
    {
        residual = function(guess) - value;
        guess = guess - residual/(function(guess + 5e-4) - function(guess - 5e-4))*1e-3;
    }
    guess
}

pub fn factorial(num: u128) -> u128
{
    (1..=num).product()
}

pub fn sinhc<T>(x: &T) -> T
where
    T: Math<T>
{
    Math::<T>::sinhc(x)
}

pub fn ln<T>(x: &T) -> T
where
    T: Math<T>
{
    Math::<T>::ln_fun(x)
}

pub fn ln_sinhc<T>(x: &T) -> T
where
    T: Math<T>
{
    Math::<T>::ln_sinhc(x)
}

pub fn langevin<T>(x: &T) -> T
where
    T: Math<T>
{
    Math::<T>::langevin(x)
}

pub fn approximate_inverse_langevin<T>(x: &T) -> T
where
    T: Math<T>
{
    Math::<T>::approximate_inverse_langevin(x)
}

pub fn inverse_langevin(y: &f64, tol: f64) -> f64
{
    let mut x = approximate_inverse_langevin(y);
    let mut residual_rel: f64 = 1.0;
    while residual_rel.abs() > tol
    {
        x = 1.0/(1.0/x.tanh() - y);
        residual_rel = langevin(&x)/y - 1.0;
    }
    return x;
}

pub trait Math<T>
{
    fn sinhc(&self) -> T;
    fn ln_fun(&self) -> T;
    fn ln_sinhc(&self) -> T;
    fn langevin(&self) -> T;
    fn approximate_inverse_langevin(&self) -> T;
}

impl Math<f64> for f64
{
    fn sinhc(&self) -> Self
    {
        self.sinh()/self
    }
    fn ln_fun(&self) -> Self
    {
        self.ln()
    }
    fn ln_sinhc(&self) -> Self
    {
        (self.sinh()/self).ln()
    }
    fn langevin(&self) -> Self
    {
        1.0/self.tanh() - 1.0/self
    }
    fn approximate_inverse_langevin(&self) -> Self
    {
        (2.14234*self.powf(3.0) - 4.22785*self.powf(2.0) + 3.0*self)/(1.0 - self)/(0.71716*self.powf(3.0) - 0.41103*self.powf(2.0) - 0.39165*self + 1.0)
    }
}
