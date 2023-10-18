extern crate polylog;
extern crate num;
use num::complex::Complex;
use num::Float;
use polylog::Li5;
mod common;


#[test]
fn special_values() {
    use num::Zero;
    let eps = 1e-15;
    let z5  = 1.0369277551433699;
    let zero = Complex::zero();

    assert_eq_complex!(zero.li5(), zero, eps);
    assert_eq_complex!(Complex::<f64>::new(1., 0.).li5(),
                       Complex::<f64>::new(z5, 0.), eps);
    assert_eq_complex!(Complex::<f64>::new(-1., 0.).li5(),
                       Complex::<f64>::new(-15./16.*z5, 0.), eps);
    assert_eq_complex!(Complex::<f64>::new(0.5, 0.).li5(),
                       Complex::<f64>::new(0.5084005792422687, 0.), eps);

    // test value that causes overflow if squared
    assert!(!Complex::new(1e300, 1.0).li5().is_infinite());
    assert!(!Complex::new(1.0, 1e300).li5().is_infinite());
    assert_eq_complex!(Complex::new(1e300, 1.0).li5(), Complex::new(-1.3105197831948743e12, 2.980481322754618e10), eps);
    assert_eq_complex!(Complex::new(1.0, 1e300).li5(), Complex::new(-1.31072310968392418e12, 1.490286896860219e10), eps);
}


#[test]
fn test_values() {
    let eps = 1e-14;
    let values = common::read_data_file("Li5.txt").unwrap();

    for &(v, li5) in values.iter() {
        assert_eq_complex!(v.li5(), li5, eps);
    }
}
