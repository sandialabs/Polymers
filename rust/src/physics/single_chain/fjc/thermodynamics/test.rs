#![cfg(test)]

mod init {

    use rand::prelude::*;
    use crate::physics::single_chain::fjc::thermodynamics::Thermodynamics;

    #[test]
    fn number_of_links()
    {
        for _ in 0..8
        {
            let number_of_links: u16 = rand::thread_rng().gen_range(8..88);
            assert_eq!(number_of_links, Thermodynamics::init(number_of_links).number_of_links);
        }
    }

    #[test]
    fn isometric()
    {
        let _ = Thermodynamics::init(8).isometric;
    }

    #[test]
    fn isotensional()
    {
        let _ = Thermodynamics::init(8).isotensional;
    }
}