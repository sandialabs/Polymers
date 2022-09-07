#![cfg(test)]

mod temporary
{
    use crate::physics::single_chain::fjc::thermodynamics::isotensional::
    {
        nondimensional_end_to_end_length_per_link,
        nondimensional_relative_gibbs_free_energy_per_link
    };

    #[test]
    fn temporary()
    {
        nondimensional_end_to_end_length_per_link(1e-3);
        nondimensional_relative_gibbs_free_energy_per_link(1e-3);
    }
}