#[cfg(feature = "python")]
pub mod py;

mod test;

/// The ideal single-chain model.
pub mod ideal;

/// The freely-jointed chain (FJC) single-chain model.
pub mod fjc;

/// The extensible freely-jointed chain (EFJC) single-chain model.
pub mod efjc;

/// The square-well freely-jointed chain (SWFJC) single-chain model.
pub mod swfjc;

/// The arbitrary link potential freely-jointed chain (uFJC) single-chain model.
pub mod ufjc;

/// The freely-rotating chain (FRC) single-chain model.
pub mod frc;

/// The extensible freely-rotating chain (EFRC) single-chain model.
pub mod efrc;

/// The worm-like chain (WLC) single-chain model.
pub mod wlc;

static ONE: f64 = 1.0;
static ZERO: f64 = 1e-6;
static POINTS: u128 = 64;
