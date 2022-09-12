#![cfg(test)]

use super::*;
use crate::physics::single_chain::test;

test::base!(Isotensional);
test::isotensional!(Isotensional);