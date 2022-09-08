#![cfg(test)]

mod init {

    use rand::prelude::*;
    use crate::physics::single_chain::fjc::FJC;

    #[test]
    fn number_of_links()
    {
        let mut rng = rand::thread_rng();
        for _ in 0..8
        {
            let number_of_links: u16 = rng.gen_range(8..88);
            assert_eq!(number_of_links, FJC::init(number_of_links).number_of_links);
        }
    }

    #[test]
    fn mechanics()
    {
        let _ = FJC::init(8).mechanics;
    }

    #[test]
    fn thermodynamics()
    {
        let _ = FJC::init(8).thermodynamics;
    }
}