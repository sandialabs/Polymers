#![cfg(test)]

use crate::physics::single_chain::test;
use crate::physics::single_chain::fjc::thermodynamics::isometric::Isometric;

test::base!(Isometric);
test::isometric!(Isometric);
