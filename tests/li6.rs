extern crate polylog;
extern crate num;
use num::complex::Complex;
use polylog::Li6;
mod common;


#[test]
fn special_values() {
    use num::Zero;
    let eps = 1e-15;
    let z6  = 1.017343061984449;
    let zero = Complex::zero();

    assert_eq_complex!(zero.li6(), zero, eps);
    assert_eq_complex!(Complex::new(1., 0.).li6(),
                       Complex::new(z6, 0.), eps);
    assert_eq_complex!(Complex::new(-1., 0.).li6(),
                       Complex::new(-31./32.*z6, 0.), eps);
    assert_eq_complex!(Complex::new(0.5, 0.).li6(),
                       Complex::new(0.5040953978039886, 0.), eps);
}


#[test]
fn test_values() {
    let eps = 1e-14;
    let values = common::read_data_file("Li6.txt").unwrap();

    for &(v, li6) in values.iter() {
        assert_eq_complex!(v.li6(), li6, eps);
    }
}
