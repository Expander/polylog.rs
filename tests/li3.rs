extern crate polylog;
extern crate num;
use num::complex::Complex;
use num::Float;
use polylog::Li3;
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

fn id1(z: Complex<f64>) -> Complex<f64> {
    z.li3() + (-z).li3() - 0.25*(z*z).li3()
}

fn id2(z: Complex<f64>) -> Complex<f64> {
    if z.norm() < std::f64::EPSILON || (z.re > 0. && z.re < 1.) {
        Complex::new(0.,0.)
    } else {
        let pi = std::f64::consts::PI;
        z.li3() - (1./z).li3() + (-z).cln().powf(3.)/6. + pi*pi/6.*(-z).cln()
    }
}

fn id3(z: Complex<f64>) -> Complex<f64> {
    if (1.0 - z).re.abs() < std::f64::EPSILON || (z.re <= 0. && z.im == 0.) {
        Complex::new(0.,0.)
    } else {
        let pi = std::f64::consts::PI;
        let z3 = 1.202056903159594;

        z.li3() + (1.-z).li3() + (1.-1./z).li3()
            - (z3 + z.cln().powf(3.)/6. + pi*pi/6.*z.cln() - 0.5*z.cln().powf(2.)*(1.-z).cln())
    }
}

#[test]
fn special_values() {
    use num::Zero;
    let pi  = std::f64::consts::PI;
    let pi2 = pi*pi;
    let eps = 1e-15;
    let ln2 = 2.0.ln();
    let z3  = 1.202056903159594;
    let phi = 0.5*((5.0).sqrt() + 1.0); // golden ratio
    let zero = Complex::zero();

    assert_eq_complex(zero.li3(), zero, eps);
    assert_eq_complex(Complex::new(1., 0.).li3(),
                      Complex::new(z3, 0.), eps);
    assert_eq_complex(Complex::new(-1., 0.).li3(),
                      Complex::new(-3./4.*z3, 0.), eps);
    assert_eq_complex(Complex::new(0.5, 0.).li3(),
                      Complex::new(ln2.powf(3.)/6. - pi2/12.*ln2 + 7./8.*z3, 0.), eps);
    assert_eq_complex(Complex::new(1./(phi*phi), 0.).li3(),
                      Complex::new(4./5.*z3 + 2./3.*phi.ln().powf(3.) - 2./15.*pi2*phi.ln(), 0.), eps);
}

#[test]
fn test_values() {
    let eps = 1e-14;
    let values = common::read_data_file("Li3.txt").unwrap();

    for &(v, li3) in values.iter() {
        assert_eq_complex(v.li3(), li3, eps);
    }
}

#[test]
fn identities() {
    use num::Zero;
    let eps = 1e-14;
    let zero = Complex::zero();
    let values = common::read_data_file("Li3.txt").unwrap();

    for &(v1, v2) in &values {
        assert_eq_complex(id1(v1), zero, eps);
        assert_eq_complex(id1(v2), zero, eps);
        assert_eq_complex(id2(v1), zero, eps);
        assert_eq_complex(id2(v2), zero, eps);
        assert_eq_complex(id3(v1), zero, eps);
        assert_eq_complex(id3(v2), zero, eps);
    }
}
