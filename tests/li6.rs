use num::complex::Complex;
use polylog::Li6;
mod common;


#[test]
fn special_values() {
    use num::Zero;
    let eps = 1e-15_f64;
    let z6  = 1.017343061984449_f64;
    let zero = Complex::zero();

    assert_eq_complex!(zero.li6(), zero, eps);
    assert_eq_complex!(Complex::<f64>::new(1.0_f64, 0.0_f64).li6(),
                       Complex::<f64>::new(z6, 0.0_f64), eps);
    assert_eq_complex!(Complex::<f64>::new(-1.0_f64, 0.0_f64).li6(),
                       Complex::<f64>::new(-31.0_f64/32.0_f64*z6, 0.0_f64), eps);
    assert_eq_complex!(Complex::<f64>::new(0.5_f64, 0.0_f64).li6(),
                       Complex::<f64>::new(0.5040953978039886_f64, 0.0_f64), eps);

    // test value that causes overflow if squared
    assert!(!Complex::new(1e300_f64, 1.0_f64).li6().is_infinite());
    assert!(!Complex::new(1.0_f64, 1e300_f64).li6().is_infinite());
    assert_eq_complex!(Complex::new(1e300_f64, 1.0_f64).li6(), Complex::new(-1.5086876165613597e14_f64, 4.11768711823317e12_f64), 2.0_f64*eps);
    assert_eq_complex!(Complex::new(1.0_f64, 1e300_f64).li6(), Complex::new(-1.5090387516918862e14_f64, 2.0589500211678e12_f64), 2.0_f64*eps);
}


#[test]
fn test_values() {
    let eps = 1e-14_f64;
    let values = common::read_data_file("Li6.txt").unwrap();

    for &(v, li6) in values.iter() {
        assert_eq_complex!(v.li6(), li6, eps);
    }
}
