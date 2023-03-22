#[cfg(feature = "python")]
pub mod py;

mod test;

// /// The uFJC single-chain model with the Lennard-Jones link potential.
// pub mod lennard_jones;

// /// The uFJC single-chain model with the log-squared link potential.
// pub mod log_squared;

/// The uFJC single-chain model with the Morse link potential.
pub mod morse;
