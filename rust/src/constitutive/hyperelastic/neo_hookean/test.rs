#![cfg(test)]

// figure out how to write this once and apply it to all hyperelastic models

mod cauchy_stress {

    use crate::constitutive::ZERO_TENSOR;
    use crate::constitutive::IDENTITY_TENSOR;
    use crate::constitutive::hyperelastic::neo_hookean::NeoHookean;

    #[test]
    fn zero() {
        let model = NeoHookean::init();
        assert_eq!(model.cauchy_stress(&*IDENTITY_TENSOR), *ZERO_TENSOR);
    }

}
