extern crate polylog;
extern crate num;
use num::complex::Complex;
use num::Float;
use polylog::Li1;
mod common;


#[test]
fn test_values() {
    let eps = 1e-14;
    let values = common::read_data_file("Li1.txt").unwrap();

    for &(v, li1) in values.iter() {
        assert_eq_complex!(v.li1(), li1, eps);

        if v.im == 0.0 {
            assert_eq_float!(v.re.li1(), li1.re, eps);
        }
    }

    assert!(1.0.li1().is_infinite());
    assert!(Complex::new(1.0, 0.0).li1().is_infinite());

    // test value that causes overflow if squared
    assert!(!Complex::new(1e300, 1.0).li1().is_infinite());
    assert!(!Complex::new(1.0, 1e300).li1().is_infinite());
    assert_eq_float!(Complex::new(1e300, 1.0).li1().re, Complex::new(-690.77552789821371, -3.14159265358979).re, eps);
    assert_eq_float!(Complex::new(1.0, 1e300).li1().re, Complex::new(-690.77552789821371, 1.5707963267948966).re, eps);
}
