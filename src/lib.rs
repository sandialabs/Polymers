#![doc = include_str!("../README.md")]

#[cfg(feature = "python")]
pub mod py;

/// Constitutive models for polymers.
pub mod constitutive;

/// Models for polymer physics.
pub mod physics;
