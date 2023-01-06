#[cfg(feature = "python")]
pub mod py;

/// Single-chain models for polymer physics.
pub mod single_chain;

/// The Boltzmann constant in units of J/(mol⋅K).
pub static BOLTZMANN_CONSTANT: f64 = 8.314462618;

/// The Planck constant in units of J⋅ns/mol.
pub static PLANCK_CONSTANT: f64 = 0.06350779923502961;