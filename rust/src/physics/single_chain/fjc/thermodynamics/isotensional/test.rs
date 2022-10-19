#![cfg(test)]

use super::*;
use crate::physics::single_chain::test;

test::base!(FJC);
test::isotensional!(FJC);

mod legendre
{
    use super::*;
    use crate::physics::single_chain::test;
    test::base!(FJCLegendre);
    // test::isotensional!(FJCLegendre);
    test::isotensionalLegendre!(FJCLegendre);
}
