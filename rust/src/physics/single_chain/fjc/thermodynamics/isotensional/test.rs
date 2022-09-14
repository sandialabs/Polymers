#![cfg(test)]

use super::*;
use crate::physics::single_chain::test_macros;

test_macros::base!(FJC);
test_macros::isotensional!(FJC);