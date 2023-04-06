#![doc = include_str!("../README.md")]

#[cfg(feature = "python")]
pub mod py;

/// Mathematical methods.
pub mod math;

/// Models for polymer physics.
pub mod physics;
