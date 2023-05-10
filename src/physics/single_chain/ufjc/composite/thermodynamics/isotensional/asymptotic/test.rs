#![cfg(test)]
use super::*;
use crate::physics::single_chain::test::Parameters;
mod base
{
    use super::*;
    use rand::Rng;
    #[test]
    fn init()
    {
        let parameters = Parameters::default();
        let _ = CUFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.number_of_bonds_minimum, parameters.bond_stiffness_reference, parameters.bond_energy_reference, parameters.bond_scission_energy_reference, parameters.bond_attempt_frequency_reference);
    }
    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            assert_eq!(number_of_links, CUFJC::init(number_of_links, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.number_of_bonds_minimum, parameters.bond_stiffness_reference, parameters.bond_energy_reference, parameters.bond_scission_energy_reference, parameters.bond_attempt_frequency_reference).number_of_links);
        }
    }
    #[test]
    fn link_length()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let link_length = parameters.link_length_reference + parameters.link_length_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(link_length, CUFJC::init(parameters.number_of_links_minimum, link_length, parameters.hinge_mass_reference, parameters.number_of_bonds_minimum, parameters.bond_stiffness_reference, parameters.bond_energy_reference, parameters.bond_scission_energy_reference, parameters.bond_attempt_frequency_reference).link_length);
        }
    }
    #[test]
    fn hinge_mass()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(hinge_mass, CUFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, hinge_mass, parameters.number_of_bonds_minimum, parameters.bond_stiffness_reference, parameters.bond_energy_reference, parameters.bond_scission_energy_reference, parameters.bond_attempt_frequency_reference).hinge_mass);
        }
    }
    #[test]
    fn number_of_bonds()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_bonds: u8 = rng.gen_range(parameters.number_of_bonds_minimum..parameters.number_of_bonds_maximum);
            assert_eq!(number_of_bonds, CUFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, number_of_bonds, parameters.bond_stiffness_reference, parameters.bond_energy_reference, parameters.bond_scission_energy_reference, parameters.bond_attempt_frequency_reference).number_of_bonds);
        }
    }
    #[test]
    fn bond_stiffness()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let bond_stiffness = parameters.bond_stiffness_reference + parameters.bond_stiffness_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(bond_stiffness, CUFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.number_of_bonds_minimum, bond_stiffness, parameters.bond_energy_reference, parameters.bond_scission_energy_reference, parameters.bond_attempt_frequency_reference).bond_stiffness);
        }
    }
    #[test]
    fn bond_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let bond_energy = parameters.bond_energy_reference + parameters.bond_energy_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(bond_energy, CUFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.number_of_bonds_minimum, parameters.bond_stiffness_reference, bond_energy, parameters.bond_scission_energy_reference, parameters.bond_attempt_frequency_reference).bond_energy);
        }
    }
    #[test]
    fn bond_scission_energy()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let bond_scission_energy = parameters.bond_scission_energy_reference + parameters.bond_scission_energy_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(bond_scission_energy, CUFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.number_of_bonds_minimum, parameters.bond_stiffness_reference, parameters.bond_energy_reference, bond_scission_energy, parameters.bond_attempt_frequency_reference).bond_scission_energy);
        }
    }
    #[test]
    fn bond_attempt_frequency()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let bond_attempt_frequency = parameters.bond_attempt_frequency_reference + parameters.bond_attempt_frequency_scale*(0.5 - rng.gen::<f64>());
            assert_eq!(bond_attempt_frequency, CUFJC::init(parameters.number_of_links_minimum, parameters.link_length_reference, parameters.hinge_mass_reference, parameters.number_of_bonds_minimum, parameters.bond_stiffness_reference, parameters.bond_energy_reference, parameters.bond_scission_energy_reference, bond_attempt_frequency).bond_attempt_frequency);
        }
    }
    #[test]
    fn all_parameters()
    {
        let mut rng = rand::thread_rng();
        let parameters = Parameters::default();
        for _ in 0..parameters.number_of_loops
        {
            let number_of_links: u8 = rng.gen_range(parameters.number_of_links_minimum..parameters.number_of_links_maximum);
            let hinge_mass = parameters.hinge_mass_reference + parameters.hinge_mass_scale*(0.5 - rng.gen::<f64>());
            let link_length = rng.gen::<f64>();
            let number_of_bonds: u8 = rng.gen_range(parameters.number_of_bonds_minimum..parameters.number_of_bonds_maximum);
            let bond_stiffness = parameters.bond_stiffness_reference + parameters.bond_stiffness_scale*(0.5 - rng.gen::<f64>());
            let bond_energy = parameters.bond_energy_reference + parameters.bond_energy_scale*(0.5 - rng.gen::<f64>());
            let bond_scission_energy = parameters.bond_scission_energy_reference + parameters.bond_scission_energy_scale*(0.5 - rng.gen::<f64>());
            let bond_attempt_frequency = parameters.bond_attempt_frequency_reference + parameters.bond_attempt_frequency_scale*(0.5 - rng.gen::<f64>());
            let model = CUFJC::init(number_of_links, link_length, hinge_mass, number_of_bonds, bond_stiffness, bond_energy, bond_scission_energy, bond_attempt_frequency);
            assert_eq!(number_of_links, model.number_of_links);
            assert_eq!(link_length, model.link_length);
            assert_eq!(hinge_mass, model.hinge_mass);
            assert_eq!(number_of_bonds, model.number_of_bonds);
            assert_eq!(bond_stiffness, model.bond_stiffness);
            assert_eq!(bond_energy, model.bond_energy);
            assert_eq!(bond_scission_energy, model.bond_scission_energy);
            assert_eq!(bond_attempt_frequency, model.bond_attempt_frequency);
        }
    }
}
