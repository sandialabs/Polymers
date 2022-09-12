#![cfg(test)]

use crate::physics::single_chain::test;
use crate::physics::single_chain::fjc::thermodynamics::Thermodynamics;

test::base!(Thermodynamics);
test::thermodynamics!(Thermodynamics);
