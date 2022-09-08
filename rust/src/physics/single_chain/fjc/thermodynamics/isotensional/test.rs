#![cfg(test)]

mod init {

    use rand::prelude::*;
    use crate::physics::single_chain::fjc::thermodynamics::isotensional::Isotensional;

    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        for _ in 0..8
        {
            let number_of_links: u16 = rng.gen_range(8..88);
            assert_eq!(number_of_links, Isotensional::init(number_of_links).number_of_links);
        }
    }
}