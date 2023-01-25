#[cfg(feature = "python")]
pub mod py;

/// Elastic constitutive models.
pub mod elastic;

// pub mod plastic;
// pub mod viscoelastic;
// pub mod viscoplastic;

fn jacobian(deformation_gradient: &[[f64; 3]; 3]) -> f64
{
      deformation_gradient[0][0]*deformation_gradient[1][1]*deformation_gradient[2][2]
    - deformation_gradient[0][0]*deformation_gradient[1][2]*deformation_gradient[2][1]
    - deformation_gradient[0][1]*deformation_gradient[1][0]*deformation_gradient[2][2]
    + deformation_gradient[0][1]*deformation_gradient[1][2]*deformation_gradient[2][0]
    + deformation_gradient[0][2]*deformation_gradient[1][0]*deformation_gradient[2][1]
    - deformation_gradient[0][2]*deformation_gradient[1][1]*deformation_gradient[2][0]
}

fn incompressible_deformation_gradient(deformation_gradient: &[[f64; 3]; 3]) -> [[f64; 3]; 3]
{
    let j_1_3 = jacobian(deformation_gradient).powf(1.0/3.0);
    [
        [
            deformation_gradient[0][0]*deformation_gradient[0][0]/j_1_3 + deformation_gradient[0][1]*deformation_gradient[0][1]/j_1_3 + deformation_gradient[0][2]*deformation_gradient[0][2]/j_1_3,
            deformation_gradient[0][0]*deformation_gradient[1][0]/j_1_3 + deformation_gradient[0][1]*deformation_gradient[1][1]/j_1_3 + deformation_gradient[0][2]*deformation_gradient[1][2]/j_1_3,
            deformation_gradient[0][0]*deformation_gradient[2][0]/j_1_3 + deformation_gradient[0][1]*deformation_gradient[2][1]/j_1_3 + deformation_gradient[0][2]*deformation_gradient[2][2]/j_1_3
        ],[
            deformation_gradient[1][0]*deformation_gradient[0][0]/j_1_3 + deformation_gradient[1][1]*deformation_gradient[0][1]/j_1_3 + deformation_gradient[1][2]*deformation_gradient[0][2]/j_1_3,
            deformation_gradient[1][0]*deformation_gradient[1][0]/j_1_3 + deformation_gradient[1][1]*deformation_gradient[1][1]/j_1_3 + deformation_gradient[1][2]*deformation_gradient[1][2]/j_1_3,
            deformation_gradient[1][0]*deformation_gradient[2][0]/j_1_3 + deformation_gradient[1][1]*deformation_gradient[2][1]/j_1_3 + deformation_gradient[1][2]*deformation_gradient[2][2]/j_1_3
        ],[
            deformation_gradient[2][0]*deformation_gradient[0][0]/j_1_3 + deformation_gradient[2][1]*deformation_gradient[0][1]/j_1_3 + deformation_gradient[2][2]*deformation_gradient[0][2]/j_1_3,
            deformation_gradient[2][0]*deformation_gradient[1][0]/j_1_3 + deformation_gradient[2][1]*deformation_gradient[1][1]/j_1_3 + deformation_gradient[2][2]*deformation_gradient[1][2]/j_1_3,
            deformation_gradient[2][0]*deformation_gradient[2][0]/j_1_3 + deformation_gradient[2][1]*deformation_gradient[2][1]/j_1_3 + deformation_gradient[2][2]*deformation_gradient[2][2]/j_1_3
        ]
    ]
}

fn left_cauchy_green(deformation_gradient: &[[f64; 3]; 3]) -> [[f64; 3]; 3]
{
    [
        [
            deformation_gradient[0][0]*deformation_gradient[0][0] + deformation_gradient[0][1]*deformation_gradient[0][1] + deformation_gradient[0][2]*deformation_gradient[0][2],
            deformation_gradient[0][0]*deformation_gradient[1][0] + deformation_gradient[0][1]*deformation_gradient[1][1] + deformation_gradient[0][2]*deformation_gradient[1][2],
            deformation_gradient[0][0]*deformation_gradient[2][0] + deformation_gradient[0][1]*deformation_gradient[2][1] + deformation_gradient[0][2]*deformation_gradient[2][2]
        ],[
            deformation_gradient[1][0]*deformation_gradient[0][0] + deformation_gradient[1][1]*deformation_gradient[0][1] + deformation_gradient[1][2]*deformation_gradient[0][2],
            deformation_gradient[1][0]*deformation_gradient[1][0] + deformation_gradient[1][1]*deformation_gradient[1][1] + deformation_gradient[1][2]*deformation_gradient[1][2],
            deformation_gradient[1][0]*deformation_gradient[2][0] + deformation_gradient[1][1]*deformation_gradient[2][1] + deformation_gradient[1][2]*deformation_gradient[2][2]
        ],[
            deformation_gradient[2][0]*deformation_gradient[0][0] + deformation_gradient[2][1]*deformation_gradient[0][1] + deformation_gradient[2][2]*deformation_gradient[0][2],
            deformation_gradient[2][0]*deformation_gradient[1][0] + deformation_gradient[2][1]*deformation_gradient[1][1] + deformation_gradient[2][2]*deformation_gradient[1][2],
            deformation_gradient[2][0]*deformation_gradient[2][0] + deformation_gradient[2][1]*deformation_gradient[2][1] + deformation_gradient[2][2]*deformation_gradient[2][2]
        ]
    ]
}

fn incompressible_left_cauchy_green(deformation_gradient: &[[f64; 3]; 3]) -> [[f64; 3]; 3]
{
    left_cauchy_green(&incompressible_deformation_gradient(deformation_gradient))
}

fn deviatoric_incompressible_left_cauchy_green(deformation_gradient: &[[f64; 3]; 3]) -> [[f64; 3]; 3]
{
    let b_bar = incompressible_left_cauchy_green(deformation_gradient);
    let b_sph = (b_bar[0][0] + b_bar[1][1] + b_bar[2][2])/3.0;
    [
        [b_bar[0][0] - b_sph, b_bar[0][1], b_bar[0][2]],
        [b_bar[1][0], b_bar[1][1] - b_sph, b_bar[1][2]],
        [b_bar[2][0], b_bar[2][1], b_bar[2][2] - b_sph],
    ]
}
