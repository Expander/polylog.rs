use num::complex::Complex;
use polylog::Li5;
mod common;


#[test]
fn special_values() {
    use num::Zero;
    let eps = 1e-15_f64;
    let z5  = 1.0369277551433699_f64;
    let zero = Complex::zero();

    assert_eq_complex!(zero.li5(), zero, eps);
    assert_eq_complex!(Complex::<f64>::new(1.0_f64, 0.0_f64).li5(),
                       Complex::<f64>::new(z5, 0.0_f64), eps);
    assert_eq_complex!(Complex::<f64>::new(-1.0_f64, 0.0_f64).li5(),
                       Complex::<f64>::new(-15.0_f64/16.0_f64*z5, 0.0_f64), eps);
    assert_eq_complex!(Complex::<f64>::new(0.5_f64, 0.0_f64).li5(),
                       Complex::<f64>::new(0.5084005792422687_f64, 0.0_f64), eps);

    // test value that causes overflow if squared
    assert!(!Complex::new(1e300_f64, 1.0_f64).li5().is_infinite());
    assert!(!Complex::new(1.0_f64, 1e300_f64).li5().is_infinite());
    assert_eq_complex!(Complex::new(1e300_f64, 1.0_f64).li5(), Complex::new(-1.3105197831948743e12_f64, 2.980481322754618e10_f64), eps);
    assert_eq_complex!(Complex::new(1.0_f64, 1e300_f64).li5(), Complex::new(-1.31072310968392418e12_f64, 1.490286896860219e10_f64), eps);
}


#[test]
fn test_values() {
    let eps = 1e-14_f64;
    let values = common::read_data_file("Li5.txt").unwrap();

    for &(v, li5) in values.iter() {
        assert_eq_complex!(v.li5(), li5, eps);
    }
}
