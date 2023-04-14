#[cfg(feature = "python")]
pub mod py;

/// The composite uFJC (CuFJC) single-chain model.
pub mod composite;

/// The uFJC single-chain model with the Lennard-Jones link potential.
pub mod lennard_jones;

/// The uFJC single-chain model with the log-squared link potential.
pub mod log_squared;

/// The uFJC single-chain model with the Morse link potential.
pub mod morse;
