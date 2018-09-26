extern crate polylog;
extern crate num;
use num::complex::Complex;
use num::Float;
use polylog::li3::Li3;

macro_rules! assert_eq_float {
    ($a:expr, $b:expr, $eps:expr) => {
        assert!(($a - $b).abs() < $eps);
    }
}

fn assert_eq_complex(a: Complex<f64>, b: Complex<f64>, eps: f64) -> () {
    assert_eq_float!(a.re, b.re, eps);
    assert_eq_float!(a.im, b.im, eps);
}

trait CLn<T> {
    fn cln(&self) -> T;
}

impl CLn<Complex<f64>> for Complex<f64> {
    fn cln(&self) -> Complex<f64> {
        let mut zf = *self;
        if zf.re == 0. { zf.re = 0. }
        if zf.im == 0. { zf.im = 0. }
        zf.ln()
    }
}

fn id1(z: Complex<f64>) -> Complex<f64> {
    z.li3() + (-z).li3() - 0.25*(z*z).li3()
}

fn id2(z: Complex<f64>) -> Complex<f64> {
    let pi = std::f64::consts::PI;

    if z.norm() < std::f64::EPSILON || (z.re > 0. && z.re < 1.) {
        Complex::new(0.,0.)
    } else {
        z.li3() - (1./z).li3() + (-z).cln().powf(3.)/6. + pi*pi/6.*(-z).cln()
    }
}

fn id3(z: Complex<f64>) -> Complex<f64> {
    let pi = std::f64::consts::PI;
    let z3 = 1.202056903159594;

    if (1.0 - z).re.abs() < 1e-10 || (z.re <= 0. && z.im == 0.) {
        Complex::new(0.,0.)
    } else {
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
    let z3  = 1.2020569031595942853997381615114;
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

/// @TODO: add Mma values
#[test]
fn test_mma_values() {
    let eps = 1e-14;
    let mma_values : Vec<(Complex<f64>, Complex<f64>)> = vec![];

    for &(v, li3) in mma_values.iter() {
        assert_eq_complex(v.li3(), li3, eps);
    }
}

/// @TODO: add Mma values
#[test]
fn test_mma_values_close_to_unity() {
    let eps = 1e-14;
    let mma_values_close_to_unity : Vec<(Complex<f64>, Complex<f64>)> = vec![];

    for &(v, li3) in mma_values_close_to_unity.iter() {
        assert_eq_complex(v.li3(), li3, eps);
    }
}

/// @TODO: add values
#[test]
fn identities() {
    use num::Zero;
    let eps = 1e-14;
    let zero = Complex::zero();
    let values : Vec<Complex<f64>> = vec![];

    for v in &values {
        assert_eq_complex(id1(*v), zero, eps);
        assert_eq_complex(id2(*v), zero, eps);
        assert_eq_complex(id3(*v), zero, eps);
    }
}
