extern crate polylog;
extern crate num;
use num::complex::Complex;
use num::Float;
use polylog::Li2;
mod common;
use common::CLn;


fn id1(z: Complex<f64>) -> Complex<f64> {
    z.li2() + (-z).li2() - 0.5*(z*z).li2()
}


fn id2(z: Complex<f64>) -> Complex<f64> {
    if z.norm() < std::f64::EPSILON || z.re < 0. {
        Complex::new(0.,0.)
    } else {
        (1.-z).li2() + (1.-1./z).li2() + 0.5*z.cln().powi(2)
    }
}


fn id3(z: Complex<f64>) -> Complex<f64> {
    let pi = std::f64::consts::PI;

    if z.norm() < 1e-10 || (z.re - 1.).abs() < 1e-10 {
        Complex::new(0.,0.)
    } else {
        z.li2() + (1.-z).li2() - pi.powi(2)/6. + z.cln()*(1.-z).cln()
    }
}


fn id4(z: Complex<f64>) -> Complex<f64> {
    if z.norm() < 1e-10 || (z.re + 1.).abs() < 1e-10 || z.re < 0. || z.im < 0. {
        Complex::new(0.,0.)
    } else {
        let pi = std::f64::consts::PI;
        (-z).li2() - (1.-z).li2() + 0.5*(1.-z*z).li2() + pi.powi(2)/12. + z.cln()*(1.+z).cln()
    }
}


fn id5(z: Complex<f64>) -> Complex<f64> {
    if z.norm() < 1e-10 || (z.re > 0. && z.re < 1.) {
        Complex::new(0.,0.)
    } else {
        let pi = std::f64::consts::PI;
        z.li2() + (1./z).li2() + pi.powi(2)/6. + 0.5*(-z).cln().powi(2)
    }
}


#[test]
fn special_values() {
    let pi = std::f64::consts::PI;
    let eps = 1e-15;

    assert_eq_float!((-1.).li2(), -pi.powi(2)/12., eps);

    assert_eq_float!((0.).li2(), 0., eps);

    assert_eq_float!((0.5).li2(), pi.powi(2)/12. - 0.5*(2.0.ln()).powi(2), eps);

    assert_eq_float!((1.).li2(), pi.powi(2)/6., eps);

    assert_eq_complex!(Complex::new(2.,0.).li2(),
                       Complex::new(pi.powi(2)/4.,0.) - Complex::i()*pi*(2.).ln(), eps);

    assert_eq_float!((-((5.).sqrt()-1.)/2.).li2(),
                     -pi.powi(2)/15. + 0.5*((((5.).sqrt()-1.)/2.).ln()).powi(2), eps);

    assert_eq_float!((-((5.).sqrt()+1.)/2.).li2(),
                     -pi.powi(2)/10. - ((((5.).sqrt()+1.)/2.).ln()).powi(2), eps);

    assert_eq_float!(((3.-(5.).sqrt())/2.).li2(),
                     pi.powi(2)/15. - ((((5.).sqrt()-1.)/2.).ln()).powi(2), eps);

    assert_eq_float!((((5.).sqrt()-1.)/2.).li2(),
                     pi.powi(2)/10. - ((((5.).sqrt()-1.)/2.).ln()).powi(2), eps);

    {
        let z = Complex::new(-1.08371e-08, 1.32716e-24);
        let li2 = z.li2();
        let li2_expected = Complex::new(-1.08370999706393160389154078878181e-8, 1.3271599928087172e-24);
        assert_eq!(li2.re, li2_expected.re);
        assert_eq!(li2.im, li2_expected.im);
    }

    // test value that causes overflow if squared
    assert!(!Complex::new(1e300, 1.0).li2().is_infinite());
    assert!(!Complex::new(1.0, 1e300).li2().is_infinite());
    assert_eq_complex!(Complex::new(1e300, 1.0).li2(), Complex::new(-238582.12510339421, 2170.13532372464), eps);
    assert_eq_complex!(Complex::new(1.0, 1e300).li2(), Complex::new(-238585.82620504462, 1085.06766186232), eps);
}


#[test]
fn special_value_identities() {
    let pi = std::f64::consts::PI;
    let eps = 1e-14;

    assert_eq_float!((1./3.).li2() - (1./9.).li2()/6.,
                     pi.powi(2)/18. - ((3.).ln()).powi(2)/6., eps);

    assert_eq_float!((-0.5).li2() + (1./9.).li2()/6.,
                     -pi.powi(2)/18. + (2.).ln()*(3.).ln()
                     - ((2.).ln()).powi(2)/2. - ((3.).ln()).powi(2)/3., eps);

    assert_eq_float!((0.25).li2() + (1./9.).li2()/3.,
                     pi.powi(2)/18. + 2.*(2.).ln()*(3.).ln()
                     - 2.*((2.).ln()).powi(2) - 2.*((3.).ln()).powi(2)/3., eps);

    assert_eq_float!((-1./3.).li2() - (1./9.).li2()/3.,
                     -pi.powi(2)/18. + ((3.).ln()).powi(2)/6., eps);

    assert_eq_float!((-1./8.).li2() + (1./9.).li2(),
                     - 0.5*((9./8.).ln()).powi(2), eps);

    assert_eq_float!(36.*(0.5).li2() - 36.*(0.25).li2()
                     - 12.*(1./8.).li2() + 6.*(1./64.).li2(),
                     pi.powi(2), 2.0*eps);
}


#[test]
fn identities() {
    use num::Zero;

    let eps = 1e-14;
    let zero = Complex::<f64>::zero();
    let omega = Complex::new(0.5, (3.).sqrt()/2.);
    let values = [
        Complex::new(0.,0.),
        Complex::new(0.5,0.),
        Complex::new(1.,0.),
        Complex::new(1.5,0.),
        Complex::new(-0.,0.),
        Complex::new(-0.5,0.),
        Complex::new(-1.,0.),
        Complex::new(-1.5,0.),
        Complex::new(-((5.).sqrt() - 1.)/2.,0.),
        Complex::new(-((5.).sqrt() + 1.)/2.,0.),
        Complex::new(((5.).sqrt() + 1.)/2.,0.),
        Complex::new(((5.).sqrt() + 3.)/2.,0.),
        omega,
        omega.powi(2),
        1. + omega,
        1./(1. + omega),
    ];

    for v in &values {
        assert_eq_complex!(id1(*v), zero, eps);
        assert_eq_complex!(id2(*v), zero, eps);
        assert_eq_complex!(id3(*v), zero, eps);
        assert_eq_complex!(id4(*v), zero, eps);
        assert_eq_complex!(id5(*v), zero, eps);
    }
}


#[test]
fn test_values() {
    let eps = 1e-14;
    let values = common::read_data_file("Li2.txt").unwrap();

    for &(v, li2) in values.iter() {
        assert_eq_complex!(v.li2(), li2, eps);

        if v.im == 0.0 {
            assert_eq_float!(v.re.li2(), li2.re, eps);
        }
    }
}
