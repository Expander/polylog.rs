extern crate polylog;
extern crate num;
use num::complex::Complex;
use polylog::Li5;
mod common;
use common::assert_eq_complex;


trait CLn<T> {
    fn cln(&self) -> T;
}

impl CLn<Complex<f64>> for Complex<f64> {
    fn cln(&self) -> Complex<f64> {
        Complex::new(
            if self.re == 0. { 0. } else { self.re },
            if self.im == 0. { 0. } else { self.im },
        ).ln()
    }
}

#[test]
fn special_values() {
    use num::Zero;
    let eps = 1e-15;
    let z5  = 1.0369277551433699;
    let zero = Complex::zero();

    assert_eq_complex(zero.li5(), zero, eps);
    assert_eq_complex(Complex::new(1., 0.).li5(),
                      Complex::new(z5, 0.), eps);
    assert_eq_complex(Complex::new(-1., 0.).li5(),
                      Complex::new(-15./16.*z5, 0.), eps);
    assert_eq_complex(Complex::new(0.5, 0.).li5(),
                      Complex::new(0.5084005792422687, 0.), eps);
}

#[test]
fn test_values() {
    let eps = 1e-14;
    let values = common::read_data_file("Li5.txt").unwrap();

    for &(v, li5) in values.iter() {
        assert_eq_complex(v.li5(), li5, eps);
    }
}
