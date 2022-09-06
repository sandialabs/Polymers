use ndarray::Array2;

pub mod test;
pub mod fjc;

pub struct SingleChain
{
}

impl SingleChain
{
    pub fn init() -> SingleChain
    {
        SingleChain
        {
        }
    }
}

pub fn sinhc<T>(nondimensional_energy: &T) -> T
where
    T: Math<T>
{
    Math::<T>::sinhc(nondimensional_energy)
}

pub fn ln_sinhc<T>(nondimensional_energy: &T) -> T
where
    T: Math<T>
{
    Math::<T>::ln_sinhc(nondimensional_energy)
}

pub fn langevin<T>(nondimensional_energy: &T) -> T
where
    T: Math<T>
{
    Math::<T>::langevin(nondimensional_energy)
}

pub trait Math<T>
{
    fn sinhc(&self) -> T;
    fn ln_sinhc(&self) -> T;
    fn langevin(&self) -> T;
}

impl Math<f64> for f64
{
    fn sinhc(&self) -> Self
    {
        self.sinh()/self
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

impl Math<Array2<f64>> for Array2<f64>
{
    fn sinhc(&self) -> Self
    {
        self.to_owned().mapv_into(|v| v.sinh())/self
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

// pub trait Math<T, A, D>
//
// impl Math<f64, _, _> for f64
//
// impl<Array<A, D>, A, D> Math for Array<A, D>
