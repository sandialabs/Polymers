#[cfg(feature = "extern")]
pub mod ex;

#[cfg(feature = "python")]
pub mod py;

mod test;

/// The worm-like chain (WLC) model thermodynamics in the isometric ensemble approximated using a Legendre transformation.
pub mod legendre;

use std::f64::consts::PI;
use crate::math::
{
    bessel_i,
    integrate_1d
};
use crate::physics::
{
    PLANCK_CONSTANT,
    BOLTZMANN_CONSTANT
};
use crate::physics::single_chain::
{
    ONE,
    ZERO,
    POINTS
};

/// The structure of the thermodynamics of the WLC model in the isometric ensemble.
pub struct WLC
{
    /// The mass of each hinge in the chain in units of kg/mol.
    pub hinge_mass: f64,

    /// The length of each link in the chain in units of nm.
    pub link_length: f64,

    /// The number of links in the chain.
    pub number_of_links: u8,

    /// The persistance length of the chain in units of nm.
    pub persistance_length: f64,

    nondimensional_persistance_length: f64,

    normalization_nondimensional_equilibrium_distribution: f64,

    /// The thermodynamic functions of the model in the isometric ensemble approximated using a Legendre transformation.
    pub legendre: self::legendre::WLC
}

/// The expected force as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, and persistance length.
pub fn force(number_of_links: &u8, link_length: &f64, persistance_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_force(number_of_links, &(persistance_length/contour_length), &(end_to_end_length/contour_length))*BOLTZMANN_CONSTANT*temperature/link_length
}

/// The expected nondimensional force as a function of the applied nondimensional end-to-end length per link, parameterized by the number of links and nondimensional persistance length.
pub fn nondimensional_force(number_of_links: &u8, nondimensional_persistance_length: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    let g2 = nondimensional_end_to_end_length_per_link.powi(2);
    let a: f64 = 14.054;
    let b: f64 = 0.473;
    let c = vec![
        vec![-0.75, 0.359375, -0.109375],
        vec![-0.5, 1.0625, -0.5625]
    ];
    let c0: f64 = 1.0 - (1.0 + (0.38/nondimensional_persistance_length.powf(0.95)).powi(-5)).powf(-0.2);
    let d: f64 = if nondimensional_persistance_length < &0.125
    {
        1.0
    }
    else
    {
        1.0 - 1.0/(0.177/(nondimensional_persistance_length - 0.111) + 6.40*(nondimensional_persistance_length - 0.111).powf(0.783))
    };
    let f = (1.0 - c0*g2)/(1.0 - g2);
    let h = nondimensional_end_to_end_length_per_link/(1.0 - b.powi(2)*g2);
    let arg = -d*nondimensional_persistance_length*a*(1.0 + b)*h;
    let sum = (0..2).collect::<Vec<usize>>().iter().map(|i| (1..4).collect::<Vec<usize>>().iter().map(|j| c[*i][j - 1]*(nondimensional_persistance_length.powi(*i as i32 - 1)*g2.powi((*j).try_into().unwrap()))).sum::<f64>()).sum::<f64>();
    let d_sum_dg = (0..2).collect::<Vec<usize>>().iter().map(|i| (1..4).collect::<Vec<usize>>().iter().map(|j| (*j as f64)*c[*i][j - 1]*(nondimensional_persistance_length.powi(*i as i32 - 1)*g2.powi((*j).try_into().unwrap()))).sum::<f64>()).sum::<f64>()*2.0/nondimensional_end_to_end_length_per_link;
    -((5.0*nondimensional_end_to_end_length_per_link*(1.0 - c0/f) + d_sum_dg + 2.0*nondimensional_end_to_end_length_per_link*sum/(1.0 - g2))/(1.0 - g2) + 2.0*b*arg *(1.0 + nondimensional_end_to_end_length_per_link*b.powi(2)*h) + arg*bessel_i(&1, &arg)/bessel_i(&0, &arg)*(1.0/nondimensional_end_to_end_length_per_link + 2.0*b.powi(2)*h))/(*number_of_links as f64)
}

/// The Helmholtz free energy as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, hinge mass, and persistance length.
pub fn helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, persistance_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, &(persistance_length/contour_length), &(end_to_end_length/contour_length), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The Helmholtz free energy per link as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, hinge mass, and persistance length.
pub fn helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, persistance_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_helmholtz_free_energy_per_link(number_of_links, link_length, hinge_mass, &(persistance_length/contour_length), &(end_to_end_length/contour_length), temperature)*BOLTZMANN_CONSTANT*temperature
}

/// The relative Helmholtz free energy as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, and persistance length.
pub fn relative_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, persistance_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_relative_helmholtz_free_energy(&(persistance_length/contour_length), &(end_to_end_length/contour_length))*BOLTZMANN_CONSTANT*temperature
}

/// The relative Helmholtz free energy per link as a function of the applied end-to-end length and temperature, parameterized by the number of links, link length, and persistance length.
pub fn relative_helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, persistance_length: &f64, end_to_end_length: &f64, temperature: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_relative_helmholtz_free_energy_per_link(number_of_links, &(persistance_length/contour_length), &(end_to_end_length/contour_length))*BOLTZMANN_CONSTANT*temperature
}

/// The nondimensional Helmholtz free energy as a function of the nondimensional end-to-end length per link and temperature, parameterized by the number of links, link length, hinge mass, and nondimensional persistance length.
pub fn nondimensional_helmholtz_free_energy(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_persistance_length: &f64, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
{
    //
    // not exactly correct unless P_eq is already normalized (see the &1.0)
    //
    let contour_length = (*number_of_links as f64)*link_length;
    -(equilibrium_distribution(number_of_links, link_length, &(contour_length*nondimensional_persistance_length), &1.0, &(contour_length*nondimensional_end_to_end_length_per_link))).ln() - ((*number_of_links as f64) - 1.0)*(4.0*(-1.0/nondimensional_persistance_length).exp().acos().sin()*PI.powi(2)*hinge_mass*link_length.powi(2)*BOLTZMANN_CONSTANT*temperature/PLANCK_CONSTANT.powi(2)).ln()
}

/// The nondimensional Helmholtz free energy per link as a function of the nondimensional end-to-end length per link and temperature, parameterized by the number of links, link length, and hinge mass, and nondimensional persistance length.
pub fn nondimensional_helmholtz_free_energy_per_link(number_of_links: &u8, link_length: &f64, hinge_mass: &f64, nondimensional_persistance_length: &f64, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
{
    nondimensional_helmholtz_free_energy(number_of_links, link_length, hinge_mass, nondimensional_persistance_length, nondimensional_end_to_end_length_per_link, temperature)/(*number_of_links as f64)
}

/// The nondimensional relative Helmholtz free energy as a function of the nondimensional end-to-end length per link, parameterized by the nondimensional persistance length.
pub fn nondimensional_relative_helmholtz_free_energy(nondimensional_persistance_length: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    (nondimensional_equilibrium_distribution(nondimensional_persistance_length, &1.0, &ZERO)/nondimensional_equilibrium_distribution(nondimensional_persistance_length, &1.0, nondimensional_end_to_end_length_per_link)).ln()
}

/// The nondimensional relative Helmholtz free energy per link as a function of the nondimensional end-to-end length per link, parameterized by the number of links and nondimensional persistance length.
pub fn nondimensional_relative_helmholtz_free_energy_per_link(number_of_links: &u8, nondimensional_persistance_length: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    nondimensional_relative_helmholtz_free_energy(nondimensional_persistance_length, nondimensional_end_to_end_length_per_link)/(*number_of_links as f64)
}

/// The equilibrium probability density of end-to-end vectors as a function of the end-to-end length, parameterized by the number of links, link length, and persistance length.
pub fn equilibrium_distribution(number_of_links: &u8, link_length: &f64, persistance_length: &f64, normalization_nondimensional_equilibrium_distribution: &f64, end_to_end_length: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_equilibrium_distribution(&(persistance_length/contour_length), normalization_nondimensional_equilibrium_distribution, &(end_to_end_length/contour_length))/contour_length.powi(3)
}

/// The nondimensional equilibrium probability density of nondimensional end-to-end vectors per link as a function of the nondimensional end-to-end length per link, parameterized by the nondimensional persistance length.
pub fn nondimensional_equilibrium_distribution(nondimensional_persistance_length: &f64, normalization_nondimensional_equilibrium_distribution: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    let g2 = nondimensional_end_to_end_length_per_link.powi(2);
    let a: f64 = 14.054;
    let b: f64 = 0.473;
    let c = vec![
        vec![-0.75, 0.359375, -0.109375],
        vec![-0.5, 1.0625, -0.5625]
    ];
    let c0: f64 = 1.0 - (1.0 + (0.38/nondimensional_persistance_length.powf(0.95)).powi(-5)).powf(-0.2);
    let d: f64;
    let e: f64;
    if nondimensional_persistance_length < &0.125
    {
        d = 1.0;
        e = (0.75/PI/nondimensional_persistance_length).powf(1.5)*(1.0 - 1.25*nondimensional_persistance_length);
    }
    else
    {
        d = 1.0 - 1.0/(0.177/(nondimensional_persistance_length - 0.111) + 6.40*(nondimensional_persistance_length - 0.111).powf(0.783));
        e = 112.04*nondimensional_persistance_length.powi(2)*(0.246/nondimensional_persistance_length - a*nondimensional_persistance_length).exp();
    }
    let f = (1.0 - c0*g2)/(1.0 - g2);
    let h = nondimensional_end_to_end_length_per_link/(1.0 - b.powi(2)*g2);
    let arg = -d*nondimensional_persistance_length*a*(1.0 + b)*h;
    let sum = (0..2).collect::<Vec<usize>>().iter().map(|i| (1..4).collect::<Vec<usize>>().iter().map(|j| c[*i][j - 1]*(nondimensional_persistance_length.powi(*i as i32 - 1)*g2.powi((*j).try_into().unwrap()))).sum::<f64>()).sum::<f64>();
    e*f.powf(2.5)*(sum/(1.0 - g2) + arg*b*nondimensional_end_to_end_length_per_link).exp()*bessel_i(&0, &arg)/normalization_nondimensional_equilibrium_distribution
}

/// The equilibrium probability density of end-to-end lengths as a function of the end-to-end length, parameterized by the number of links, link length, and persistance length.
pub fn equilibrium_radial_distribution(number_of_links: &u8, link_length: &f64, persistance_length: &f64, normalization_nondimensional_equilibrium_distribution: &f64, end_to_end_length: &f64) -> f64
{
    let contour_length = (*number_of_links as f64)*link_length;
    nondimensional_equilibrium_radial_distribution(&(persistance_length/contour_length), normalization_nondimensional_equilibrium_distribution, &(end_to_end_length/contour_length))/contour_length
}

/// The nondimensional equilibrium probability density of nondimensional end-to-end lengths per link as a function of the nondimensional end-to-end length per link, parameterized by the nondimensional persistance length.
pub fn nondimensional_equilibrium_radial_distribution(nondimensional_persistance_length: &f64, normalization_nondimensional_equilibrium_distribution: &f64, nondimensional_end_to_end_length_per_link: &f64) -> f64
{
    4.0*PI*nondimensional_end_to_end_length_per_link.powi(2)*nondimensional_equilibrium_distribution(nondimensional_persistance_length, normalization_nondimensional_equilibrium_distribution, nondimensional_end_to_end_length_per_link)
}

/// The implemented functionality of the thermodynamics of the WLC model in the isometric ensemble.
impl WLC
{
    /// Initializes and returns an instance of the thermodynamics of the WLC model in the isometric ensemble.
    pub fn init(number_of_links: u8, link_length: f64, hinge_mass: f64, persistance_length: f64) -> Self
    {
        let nondimensional_persistance_length = persistance_length/(number_of_links as f64)/link_length;
        let normalization = integrate_1d(&|nondimensional_end_to_end_length_per_link: &f64| nondimensional_equilibrium_radial_distribution(&nondimensional_persistance_length, &1.0, nondimensional_end_to_end_length_per_link), &ZERO, &ONE, &POINTS);
        WLC
        {
            hinge_mass,
            link_length,
            number_of_links,
            persistance_length,
            nondimensional_persistance_length: nondimensional_persistance_length,
            normalization_nondimensional_equilibrium_distribution: normalization,
            legendre: self::legendre::WLC::init(number_of_links, link_length, hinge_mass, persistance_length),
        }
    }
    /// The expected force as a function of the applied end-to-end length and temperature.
    pub fn force(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        force(&self.number_of_links, &self.link_length, &self.persistance_length, end_to_end_length, temperature)
    }
    /// The expected nondimensional force as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_force(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_force(&self.number_of_links, &self.nondimensional_persistance_length, nondimensional_end_to_end_length_per_link)
    }
    /// The Helmholtz free energy as a function of the applied end-to-end length and temperature.
    pub fn helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.persistance_length, end_to_end_length, temperature)
    }
    /// The Helmholtz free energy per link as a function of the applied end-to-end length and temperature.
    pub fn helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.persistance_length, end_to_end_length, temperature)
    }
    /// The relative Helmholtz free energy as a function of the applied end-to-end length and temperature.
    pub fn relative_helmholtz_free_energy(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.persistance_length, end_to_end_length, temperature)
    }
    /// The relative Helmholtz free energy per link as a function of the applied end-to-end length and temperature.
    pub fn relative_helmholtz_free_energy_per_link(&self, end_to_end_length: &f64, temperature: &f64) -> f64
    {
        relative_helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.persistance_length, end_to_end_length, temperature)
    }
    /// The nondimensional Helmholtz free energy as a function of the applied nondimensional end-to-end length per link and temperature.
    pub fn nondimensional_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.nondimensional_persistance_length, nondimensional_end_to_end_length_per_link, temperature)
    }
    /// The nondimensional Helmholtz free energy per link as a function of the applied nondimensional end-to-end length per link and temperature.
    pub fn nondimensional_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64, temperature: &f64) -> f64
    {
        nondimensional_helmholtz_free_energy_per_link(&self.number_of_links, &self.link_length, &self.hinge_mass, &self.nondimensional_persistance_length, nondimensional_end_to_end_length_per_link, temperature)
    }
    /// The nondimensional relative Helmholtz free energy as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_relative_helmholtz_free_energy(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy(&self.nondimensional_persistance_length, nondimensional_end_to_end_length_per_link)
    }
    /// The nondimensional relative Helmholtz free energy per link as a function of the applied nondimensional end-to-end length per link.
    pub fn nondimensional_relative_helmholtz_free_energy_per_link(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_relative_helmholtz_free_energy_per_link(&self.number_of_links, &self.nondimensional_persistance_length, nondimensional_end_to_end_length_per_link)
    }
    /// The equilibrium probability density of end-to-end vectors as a function of the end-to-end length.
    pub fn equilibrium_distribution(&self, end_to_end_length: &f64) -> f64
    {
        equilibrium_distribution(&self.number_of_links, &self.link_length, &self.persistance_length, &self.normalization_nondimensional_equilibrium_distribution, end_to_end_length)
    }
    /// The nondimensional equilibrium probability density of nondimensional end-to-end vectors per link as a function of the nondimensional end-to-end length per link.
    pub fn nondimensional_equilibrium_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_equilibrium_distribution(&self.nondimensional_persistance_length, &self.normalization_nondimensional_equilibrium_distribution, nondimensional_end_to_end_length_per_link)
    }
    /// The equilibrium probability density of end-to-end lengths as a function of the end-to-end length.
    pub fn equilibrium_radial_distribution(&self, end_to_end_length: &f64) -> f64
    {
        equilibrium_radial_distribution(&self.number_of_links, &self.link_length, &self.persistance_length, &self.normalization_nondimensional_equilibrium_distribution, end_to_end_length)
    }
    /// The nondimensional equilibrium probability density of nondimensional end-to-end lengths per link as a function of the nondimensional end-to-end length per link.
    pub fn nondimensional_equilibrium_radial_distribution(&self, nondimensional_end_to_end_length_per_link: &f64) -> f64
    {
        nondimensional_equilibrium_radial_distribution(&self.nondimensional_persistance_length, &self.normalization_nondimensional_equilibrium_distribution, nondimensional_end_to_end_length_per_link)
    }
}
