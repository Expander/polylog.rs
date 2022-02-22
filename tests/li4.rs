extern crate polylog;
extern crate num;
use num::complex::Complex;
use polylog::Li4;
mod common;


#[test]
fn special_values() {
    use num::Zero;
    let eps = 1e-15;
    let z4  = 1.082323233711138;
    let zero = Complex::zero();

    assert_eq_complex!(zero.li4(), zero, eps);
    assert_eq_complex!(Complex::<f64>::new(1., 0.).li4(),
                       Complex::<f64>::new(z4, 0.), eps);
    assert_eq_complex!(Complex::<f64>::new(-1., 0.).li4(),
                       Complex::<f64>::new(-7./8.*z4, 0.), eps);
    assert_eq_complex!(Complex::<f64>::new(0.5, 0.).li4(),
                       Complex::<f64>::new(0.5174790616738994, 0.), eps);
}


#[test]
fn test_values() {
    let eps = 1e-14;
    let values = common::read_data_file("Li4.txt").unwrap();

    for &(v, li4) in values.iter() {
        assert_eq_complex!(v.li4(), li4, eps);

        if v.im == 0.0 {
            assert_eq_float!(v.re.li4(), li4.re, eps);
        }
    }
}
