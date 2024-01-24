use num::complex::Complex;
use polylog::Li0;
mod common;


#[test]
fn test_values() {
    let eps = 1e-14_f64;
    let values = common::read_data_file("Li0.txt").unwrap();

    for &(v, li0) in values.iter() {
        assert_eq_complex!(v.li0(), li0, eps);

        if v.im == 0.0_f64 {
            assert_eq_float!(v.re.li0(), li0.re, eps);
        }
    }

    assert!(1.0_f64.li0().is_infinite());
    assert!(Complex::new(1.0_f64, 0.0_f64).li0().is_nan());

    // test value that causes overflow if squared
    assert!(!Complex::new(1e300_f64, 1.0_f64).li0().is_infinite());
    assert!(!Complex::new(1.0_f64, 1e300_f64).li0().is_infinite());
    assert_eq_complex!(Complex::new(1e300_f64, 1.0_f64).li0(), Complex::new(-1.0_f64, 0.0_f64), eps);
    assert_eq_complex!(Complex::new(1.0_f64, 1e300_f64).li0(), Complex::new(-1.0_f64, 0.0_f64), eps);
}
