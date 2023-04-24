#[cfg(feature = "extern")]
pub mod ex;

#[cfg(feature = "python")]
pub mod py;

mod test;

/// The worm-like chain (WLC) model thermodynamics in the isometric ensemble approximated using a Legendre transformation.
pub mod legendre;

/// The structure of the thermodynamics of the WLC model in the isometric ensemble.
pub struct WLC
{
    /// The contour length of the chain in units of nm.
    pub chain_length: f64,

    /// The persistance length of the chain in units of nm.
    pub persistance_length: f64,

    /// The thermodynamic functions of the model in the isometric ensemble approximated using a Legendre transformation.
    pub legendre: self::legendre::WLC
}

/// The implemented functionality of the thermodynamics of the WLC model in the isometric ensemble.
impl WLC
{
    /// Initializes and returns an instance of the thermodynamics of the WLC model in the isometric ensemble.
    pub fn init(chain_length: f64, persistance_length: f64) -> Self
    {
        WLC
        {
            chain_length,
            persistance_length,
            legendre: self::legendre::WLC::init(chain_length, persistance_length),
        }
    }
}
