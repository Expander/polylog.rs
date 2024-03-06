use num::complex::Complex;
use polylog::Li2;
mod common;
use common::CLn;


fn id1(z: Complex<f64>) -> Complex<f64> {
    z.li2() + (-z).li2() - 0.5_f64*(z*z).li2()
}


fn id2(z: Complex<f64>) -> Complex<f64> {
    if z.norm() < std::f64::EPSILON || z.re < 0.0_f64 {
        Complex::new(0.0_f64, 0.0_f64)
    } else {
        (1.0_f64 - z).li2() + (1.0_f64 - 1.0_f64/z).li2() + 0.5_f64*z.cln().powi(2)
    }
}


fn id3(z: Complex<f64>) -> Complex<f64> {
    let pi = std::f64::consts::PI;

    if z.norm() < 1e-10_f64 || (z.re - 1.0_f64).abs() < 1e-10_f64 {
        Complex::new(0.0_f64, 0.0_f64)
    } else {
        z.li2() + (1.0_f64 - z).li2() - pi.powi(2)/6.0_f64 + z.cln()*(1.0_f64 - z).cln()
    }
}


fn id4(z: Complex<f64>) -> Complex<f64> {
    if z.norm() < 1e-10_f64 || (z.re + 1.0_f64).abs() < 1e-10_f64 || z.re < 0.0_f64 || z.im < 0.0_f64 {
        Complex::new(0.0_f64, 0.0_f64)
    } else {
        let pi = std::f64::consts::PI;
        (-z).li2() - (1.0_f64 - z).li2() + 0.5_f64*(1.0_f64 - z*z).li2() + pi.powi(2)/12.0_f64 + z.cln()*(1.0_f64 + z).cln()
    }
}


fn id5(z: Complex<f64>) -> Complex<f64> {
    if z.norm() < 1e-10_f64 || (z.re > 0.0_f64 && z.re < 1.0_f64) {
        Complex::new(0.0_f64, 0.0_f64)
    } else {
        let pi = std::f64::consts::PI;
        z.li2() + (1.0_f64/z).li2() + pi.powi(2)/6.0_f64 + 0.5_f64*(-z).cln().powi(2)
    }
}


#[test]
fn special_values() {
    let pi = std::f64::consts::PI;
    let eps = 1e-15_f64;

    assert_eq_float!((-1.0_f64).li2(), -pi.powi(2)/12.0_f64, eps);

    assert_eq_float!((0.0_f64).li2(), 0.0_f64, eps);

    assert_eq_float!((0.5_f64).li2(), pi.powi(2)/12.0_f64 - 0.5_f64*(2.0_f64.ln()).powi(2), eps);

    assert_eq_float!((1.0_f64).li2(), pi.powi(2)/6.0_f64, eps);

    assert_eq_complex!(Complex::new(2.0_f64, 0.0_f64).li2(),
                       Complex::new(pi.powi(2)/4.0_f64, 0.0_f64) - Complex::i()*pi*2.0_f64.ln(), eps);

    assert_eq_float!((-(5.0_f64.sqrt() - 1.0_f64)/2.0_f64).li2(),
                     -pi.powi(2)/15.0_f64 + 0.5_f64*(((5.0_f64.sqrt() - 1.0_f64)/2.0_f64).ln()).powi(2), eps);

    assert_eq_float!((-(5.0_f64.sqrt() + 1.0_f64)/2.0_f64).li2(),
                     -pi.powi(2)/10.0_f64 - (((5.0_f64.sqrt() + 1.0_f64)/2.0_f64).ln()).powi(2), eps);

    assert_eq_float!(((3.0_f64 - 5.0_f64.sqrt())/2.0_f64).li2(),
                     pi.powi(2)/15.0_f64 - (((5.0_f64.sqrt() - 1.0_f64)/2.0_f64).ln()).powi(2), eps);

    assert_eq_float!(((5.0_f64.sqrt() - 1.0_f64)/2.0_f64).li2(),
                     pi.powi(2)/10.0_f64 - (((5.0_f64.sqrt() - 1.0_f64)/2.0_f64).ln()).powi(2), eps);

    {
        let z = Complex::new(-1.08371e-08_f64, 1.32716e-24_f64);
        let li2 = z.li2();
        let li2_expected = Complex::new(-1.08370999706393160389154078878181e-8_f64, 1.3271599928087172e-24_f64);
        assert_eq!(li2.re, li2_expected.re);
        assert_eq!(li2.im, li2_expected.im);
    }

    // test value that causes overflow if squared
    assert!(!Complex::new(1e300_f64, 1.0_f64).li2().is_infinite());
    assert!(!Complex::new(1.0_f64, 1e300_f64).li2().is_infinite());
    assert_eq_complex!(Complex::new(1e300_f64, 1.0_f64).li2(), Complex::new(-238582.12510339421_f64, 2170.13532372464_f64), eps);
    assert_eq_complex!(Complex::new(1.0_f64, 1e300_f64).li2(), Complex::new(-238585.82620504462_f64, 1085.06766186232_f64), eps);
}


#[test]
fn special_value_identities() {
    let pi = std::f64::consts::PI;
    let eps = 1e-14_f64;

    assert_eq_float!((1.0_f64/3.0_f64).li2() - (1.0_f64/9.0_f64).li2()/6.0_f64,
                     pi.powi(2)/18.0_f64 - ((3.0_f64).ln()).powi(2)/6.0_f64, eps);

    assert_eq_float!((-0.5_f64).li2() + (1.0_f64/9.0_f64).li2()/6.0_f64,
                     -pi.powi(2)/18.0_f64 + (2.0_f64).ln()*(3.0_f64).ln()
                     - ((2.0_f64).ln()).powi(2)/2.0_f64 - ((3.0_f64).ln()).powi(2)/3.0_f64, eps);

    assert_eq_float!((0.25_f64).li2() + (1.0_f64/9.0_f64).li2()/3.0_f64,
                     pi.powi(2)/18.0_f64 + 2.0_f64*(2.0_f64).ln()*(3.0_f64).ln()
                     - 2.0_f64*((2.0_f64).ln()).powi(2) - 2.0_f64*((3.0_f64).ln()).powi(2)/3.0_f64, eps);

    assert_eq_float!((-1.0_f64/3.0_f64).li2() - (1.0_f64/9.0_f64).li2()/3.0_f64,
                     -pi.powi(2)/18.0_f64 + ((3.0_f64).ln()).powi(2)/6.0_f64, eps);

    assert_eq_float!((-1.0_f64/8.0_f64).li2() + (1.0_f64/9.0_f64).li2(),
                     - 0.5_f64*((9.0_f64/8.0_f64).ln()).powi(2), eps);

    assert_eq_float!(36.0_f64*(0.5_f64).li2() - 36.0_f64*(0.25_f64).li2()
                     - 12.0_f64*(1.0_f64/8.0_f64).li2() + 6.0_f64*(1.0_f64/64.0_f64).li2(),
                     pi.powi(2), 2.0_f64*eps);
}


#[test]
fn identities() {
    use num::Zero;

    let eps = 1e-14_f64;
    let zero = Complex::<f64>::zero();
    let omega = Complex::new(0.5_f64, (3.0_f64).sqrt()/2.0_f64);
    let values = [
        Complex::new(0.0_f64, 0.0_f64),
        Complex::new(0.5_f64, 0.0_f64),
        Complex::new(1.0_f64, 0.0_f64),
        Complex::new(1.5_f64, 0.0_f64),
        Complex::new(-0.0_f64, 0.0_f64),
        Complex::new(-0.5_f64, 0.0_f64),
        Complex::new(-1.0_f64, 0.0_f64),
        Complex::new(-1.5_f64, 0.0_f64),
        Complex::new(-((5.0_f64).sqrt() - 1.0_f64)/2.0_f64, 0.0_f64),
        Complex::new(-((5.0_f64).sqrt() + 1.0_f64)/2.0_f64, 0.0_f64),
        Complex::new(((5.0_f64).sqrt() + 1.0_f64)/2.0_f64, 0.0_f64),
        Complex::new(((5.0_f64).sqrt() + 3.0_f64)/2.0_f64, 0.0_f64),
        omega,
        omega.powi(2),
        1.0_f64 + omega,
        1.0_f64/(1.0_f64 + omega),
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
fn test_values_f32() {
    let eps = 10.0_f32*std::f32::EPSILON;
    let values = common::read_data_file::<f32>("Li2.txt").unwrap();

    for &(v, li2) in values.iter() {
        assert_eq_complex!(v.li2(), li2, eps);

        if v.im == 0.0_f32 {
            assert_eq_float!(v.re.li2(), li2.re, eps);
        }
    }
}


#[test]
fn test_values_f64() {
    let eps = 10.0_f64*std::f64::EPSILON;
    let values = common::read_data_file::<f64>("Li2.txt").unwrap();

    for &(v, li2) in values.iter() {
        assert_eq_complex!(v.li2(), li2, eps);

        if v.im == 0.0_f64 {
            assert_eq_float!(v.re.li2(), li2.re, eps);
        }
    }
}


#[test]
fn test_signed_zero() {
    let pz32 = 0.0_f32;
    let nz32 = -0.0_f32;
    let pz64 = 0.0_f64;
    let nz64 = -0.0_f64;

    assert!(pz32.li2().is_sign_positive());
    assert!(nz32.li2().is_sign_negative());

    assert!(pz64.li2().is_sign_positive());
    assert!(nz64.li2().is_sign_negative());

    assert!(Complex::new(pz32, pz32).li2().re.is_sign_positive());
    assert!(Complex::new(pz32, pz32).li2().im.is_sign_positive());
    assert!(Complex::new(pz32, nz32).li2().re.is_sign_positive());
    assert!(Complex::new(pz32, nz32).li2().im.is_sign_negative());
    assert!(Complex::new(nz32, pz32).li2().re.is_sign_negative());
    assert!(Complex::new(nz32, pz32).li2().im.is_sign_positive());
    assert!(Complex::new(nz32, nz32).li2().re.is_sign_negative());
    assert!(Complex::new(nz32, nz32).li2().im.is_sign_negative());

    assert!(Complex::new(pz64, pz64).li2().re.is_sign_positive());
    assert!(Complex::new(pz64, pz64).li2().im.is_sign_positive());
    assert!(Complex::new(pz64, nz64).li2().re.is_sign_positive());
    assert!(Complex::new(pz64, nz64).li2().im.is_sign_negative());
    assert!(Complex::new(nz64, pz64).li2().re.is_sign_negative());
    assert!(Complex::new(nz64, pz64).li2().im.is_sign_positive());
    assert!(Complex::new(nz64, nz64).li2().re.is_sign_negative());
    assert!(Complex::new(nz64, nz64).li2().im.is_sign_negative());
}
