#![cfg(test)]

mod init {

    use rand::prelude::*;
    use crate::physics::single_chain::fjc::mechanics::Mechanics;

    #[test]
    fn number_of_links()
    {
        for _ in 0..8
        {
            let number_of_links: u16 = rand::thread_rng().gen_range(8..88);
            assert_eq!(number_of_links, Mechanics::init(number_of_links).number_of_links);
        }
    }
}

mod random_configuration {

    use rand::prelude::*;
    use crate::physics::single_chain::fjc::mechanics::Mechanics;

    #[test]
    fn number_of_configurations()
    {
        for _ in 0..8
        {
            let number_of_links: u16 = rand::thread_rng().gen_range(8..88);
            let model = Mechanics::init(number_of_links);
            assert_eq!(number_of_links + 1, model.random_configuration().len() as u16);
        }
    }
}