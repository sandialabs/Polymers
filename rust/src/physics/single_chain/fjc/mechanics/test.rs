#![cfg(test)]

use crate::physics::single_chain::test_macros;
use crate::physics::single_chain::fjc::mechanics::Mechanics;

test_macros::base!(Mechanics);
test_macros::mechanics!(Mechanics);
