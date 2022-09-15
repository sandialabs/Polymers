use ndarray::{
    Array,
    Dimension
};

pub mod test;

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

pub trait Math<T>
{
    fn sinhc(&self) -> T;
    fn ln_fun(&self) -> T;
    fn ln_sinhc(&self) -> T;
    fn langevin(&self) -> T;
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
}

impl<D> Math<Array<f64, D>> for Array<f64, D>
where
    D: Dimension
{
    fn sinhc(&self) -> Self
    {
        self.to_owned().mapv_into(|v| v.sinh())/self
    }
    fn ln_fun(&self) -> Self
    {
        self.to_owned().mapv_into(|v| v.ln())
    }
    fn ln_sinhc(&self) -> Self
    {
        (self.to_owned().mapv_into(|v| v.sinh())/self).mapv_into(|v| v.ln())
    }
    fn langevin(&self) -> Self
    {
        1.0/self.to_owned().mapv_into(|v| v.tanh()) - 1.0/self
    }
}