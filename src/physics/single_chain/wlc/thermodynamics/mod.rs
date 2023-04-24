#[cfg(feature = "python")]
pub mod py;

mod test;

/// The worm-like chain (WLC) model thermodynamics in the isometric ensemble.
pub mod isometric;

/// The structure of the thermodynamics of the WLC model.
pub struct WLC
{
    /// The contour length of the chain in units of nm.
    pub chain_length: f64,

    /// The persistance length of the chain in units of nm.
    pub persistance_length: f64,

    /// The thermodynamic functions of the model in the isometric ensemble.
    pub isometric: self::isometric::WLC
}

/// The implemented functionality of the thermodynamics of the WLC model.
impl WLC
{
    /// Initializes and returns an instance of the thermodynamics of the WLC model.
    pub fn init(chain_length: f64, persistance_length: f64) -> Self
    {
        WLC
        {
            chain_length,
            persistance_length,
            isometric: self::isometric::WLC::init(chain_length, persistance_length),
        }
    }
}
