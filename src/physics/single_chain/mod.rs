#[cfg(feature = "python")]
pub mod py;

mod test;

/// The ideal single-chain model.
pub mod ideal;

/// The freely-jointed chain (FJC) single-chain model.
pub mod fjc;

/// The extensible freely-jointed chain (EFJC) single-chain model.
pub mod efjc;

/// The square-well freely-jointed chain (EFJC) single-chain model.
pub mod swfjc;

static ONE: f64 = 1.0;
static ZERO: f64 = 1e-6;
static POINTS: u128 = 100;
