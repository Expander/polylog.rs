extern crate polylog;
extern crate num;
use num::complex::Complex;
use polylog::Li0;
mod common;


#[test]
fn test_values() {
    let eps = 1e-14;
    let values = common::read_data_file("Li0.txt").unwrap();

    for &(v, li0) in values.iter() {
        assert_eq_complex!(v.li0(), li0, eps);

        if v.im == 0.0 {
            assert_eq_float!(v.re.li0(), li0.re, eps);
        }
    }

    assert!(1.0.li0().is_infinite());
    assert!(Complex::new(1.0, 0.0).li0().is_nan());
}
