use num::complex::Complex;
use polylog::Li4;
mod common;


#[test]
fn special_values() {
    use num::Zero;
    let eps = 1e-15_f64;
    let z4  = 1.082323233711138_f64;
    let zero = Complex::zero();

    assert_eq_complex!(zero.li4(), zero, eps);
    assert_eq_complex!(Complex::<f64>::new(1.0_f64, 0.0_f64).li4(),
                       Complex::<f64>::new(z4, 0.0_f64), eps);
    assert_eq_complex!(Complex::<f64>::new(-1.0_f64, 0.0_f64).li4(),
                       Complex::<f64>::new(-7.0_f64/8.0_f64*z4, 0.0_f64), eps);
    assert_eq_complex!(Complex::<f64>::new(0.5_f64, 0.0_f64).li4(),
                       Complex::<f64>::new(0.5174790616738994_f64, 0.0_f64), eps);

    // test value that causes overflow if squared
    assert!(!Complex::new(1e300_f64, 1.0_f64).li4().is_infinite());
    assert!(!Complex::new(1.0_f64, 1e300_f64).li4().is_infinite());
    assert_eq_complex!(Complex::new(1e300_f64, 1.0_f64).li4(), Complex::new(-9.4863817894708364e9_f64, 1.725875455850714e8_f64), eps);
    assert_eq_complex!(Complex::new(1.0_f64, 1e300_f64).li4(), Complex::new(-9.4872648206269765e9_f64, 8.62951114411071e7_f64), eps);
}


#[test]
fn test_values() {
    let eps = 1e-14_f64;
    let values = common::read_data_file("Li4.txt").unwrap();

    for &(v, li4) in values.iter() {
        assert_eq_complex!(v.li4(), li4, eps);

        if v.im == 0.0_f64 {
            assert_eq_float!(v.re.li4(), li4.re, eps);
        }
    }
}
