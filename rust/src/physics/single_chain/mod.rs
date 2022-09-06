use ndarray::Array2;

pub mod test;
pub mod fjc;

pub struct SingleChain {
}

impl SingleChain {
    pub fn init() -> SingleChain {
        SingleChain {
        }
    }
}

pub trait Math {
    fn sinh(&self) -> Self;
    fn langevin(&self) -> Self;
}

impl Math for f64 {
    fn sinh(&self) -> Self {
        self.sinh()
    }
    fn langevin(&self) -> Self {
        1.0/self.tanh() - 1.0/self
    }
}

impl Math for Array2<f64> {
    fn sinh(&self) -> Self {
        self.to_owned().mapv_into(|v| v.sinh())
    }
    fn langevin(&self) -> Self {
        1.0/self.to_owned().mapv_into(|v| v.tanh()) - 1.0/self
    }
}

// impl<A, D> Math for Array<A, D> {
//     fn sinh(&self) -> Self {
//         self.to_owned().mapv_into(|v| v.sinh())
//     }
//     fn langevin(&self) -> Self {
//         self.to_owned().mapv_into(|v| v.sinh())
//     }
// }
