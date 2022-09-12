#![cfg(test)]

use crate::physics::single_chain::test;
use crate::physics::single_chain::fjc::FJC;

test::base!(FJC);
test::single_chain!(FJC);
