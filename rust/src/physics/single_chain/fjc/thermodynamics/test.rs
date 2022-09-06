#![cfg(test)]

mod temporary {

    use crate::physics::single_chain::fjc::thermodynamics::{
        nondimensional_end_to_end_length_per_link,
        nondimensional_gibbs_free_energy_per_link
    };

    #[test]
    fn temporary() {
        nondimensional_end_to_end_length_per_link(0.0);
        nondimensional_gibbs_free_energy_per_link(0.0);
    }

}
