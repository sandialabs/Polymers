extern crate polymers;

use polymers::physics::single_chain::fjc::FJC;

#[test]
fn temp()
{
    let _ = FJC::init(8, 1.0, 1.0).thermodynamics.isometric.force(&0.5, &300.0);
}
