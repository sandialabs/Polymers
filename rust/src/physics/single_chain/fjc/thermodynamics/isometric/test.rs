#![cfg(test)]

use super::*;
use crate::physics::single_chain::test;

test::base!(FJC);
test::isometric!(FJC);
mod legendre
{
    use super::*;
    use crate::physics::single_chain::test;
    test::base!(FJCLegendre);
    test::isometric!(FJCLegendre);
    test::isometricLegendre!(FJCLegendre);
}
