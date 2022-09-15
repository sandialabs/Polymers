#![cfg(test)]

use crate::physics::single_chain::test;
use crate::physics::single_chain::fjc::mechanics::Mechanics;

test::base!(Mechanics);
test::mechanics!(Mechanics);
