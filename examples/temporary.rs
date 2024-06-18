fn main() {
    // gauss is way faster if can justify using it
    // seems like GL-GL reproduces past results, but seems like GL-Gauss is still off, why?
    let num: usize = 33;
    let mut lambda_11: f64;
    let model = polymers::constitutive::hyperelastic::BucheSilberstein::init(3, 50.0, 25);
    for index in 0..num {
        lambda_11 = 1.0 + 2.5 * (index as f64) / (num as f64 - 1.0);
        println!("{:?}", (lambda_11, model.uniaxial_tension(&lambda_11)))
    }
}