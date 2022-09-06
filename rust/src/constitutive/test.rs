#![cfg(test)]

mod matrix_multiply {

    use crate::constitutive::ZERO_TENSOR;
    use crate::constitutive::IDENTITY_TENSOR;
    use crate::constitutive::matrix_multiply;

    #[test]
    fn zero_tensors() {
        assert_eq!(matrix_multiply(&*ZERO_TENSOR, &*ZERO_TENSOR), *ZERO_TENSOR);
    }

    #[test]
    fn identity_tensors() {
        assert_eq!(matrix_multiply(&*IDENTITY_TENSOR, &*IDENTITY_TENSOR), *IDENTITY_TENSOR);
    }

    #[test]
    fn zero_and_identity_tensors() {
        assert_eq!(matrix_multiply(&*ZERO_TENSOR, &*IDENTITY_TENSOR), *ZERO_TENSOR);
    }

}
