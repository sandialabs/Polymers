#![cfg(test)]

use crate::physics::single_chain::test_macros;
use crate::physics::single_chain::fjc::thermodynamics::Thermodynamics;

test_macros::base!(Thermodynamics);
test_macros::thermodynamics!(Thermodynamics);
