use num::complex::Complex;
use num::Float;
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

    // test value that causes overflow if squared
    assert!(!Complex::new(1e300, 1.0).li4().is_infinite());
    assert!(!Complex::new(1.0, 1e300).li4().is_infinite());
    assert_eq_complex!(Complex::new(1e300, 1.0).li4(), Complex::new(-9.4863817894708364e9, 1.725875455850714e8), eps);
    assert_eq_complex!(Complex::new(1.0, 1e300).li4(), Complex::new(-9.4872648206269765e9, 8.62951114411071e7), eps);
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
