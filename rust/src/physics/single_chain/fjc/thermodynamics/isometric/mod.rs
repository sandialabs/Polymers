use crate::math::
{
    Math,
    approximate_inverse_langevin
};
use crate::physics::single_chain::
{
    Isometric,
    IsometricLegendre
};

pub mod test;

pub struct FJC
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u16,
    pub number_of_links_f64: f64,
    pub legendre: FJCLegendre
}

impl Isometric for FJC
{
    fn init(number_of_links: u16, link_length: f64, hinge_mass: f64) -> FJC
    {
        FJC
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_links_f64: number_of_links as f64,
            legendre: FJCLegendre::init(number_of_links, link_length, hinge_mass)
        }
    }
}

pub struct FJCLegendre
{
    pub hinge_mass: f64,
    pub link_length: f64,
    pub number_of_links: u16,
    pub number_of_links_f64: f64
}

impl IsometricLegendre for FJCLegendre
{
    fn init(number_of_links: u16, link_length: f64, hinge_mass: f64) -> FJCLegendre
    {
        FJCLegendre
        {
            hinge_mass,
            link_length,
            number_of_links,
            number_of_links_f64: number_of_links as f64
        }
    }
    fn nondimensional_force<T>(&self, nondimensional_end_to_end_length_per_link: &T) -> T
    where T:
        Math<T> +
        std::marker::Copy +
        std::ops::Neg<Output = T> +
        std::ops::Mul<T, Output = T> +
        std::ops::Sub<T, Output = T> +
        std::ops::Add<f64, Output = T> +
        std::ops::Div<f64, Output = T> +
        std::ops::Mul<f64, Output = T>
    {
        approximate_inverse_langevin(nondimensional_end_to_end_length_per_link)
    }
}
