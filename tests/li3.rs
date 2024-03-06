use num::complex::Complex;
use polylog::Li3;
mod common;
use common::CLn;


fn id1(z: Complex<f64>) -> Complex<f64> {
    z.li3() + (-z).li3() - 0.25_f64*(z*z).li3()
}


fn id2(z: Complex<f64>) -> Complex<f64> {
    if z.norm() < std::f64::EPSILON || (z.re > 0.0_f64 && z.re < 1.0_f64) {
        Complex::new(0.0_f64, 0.0_f64)
    } else {
        let pi = std::f64::consts::PI;
        z.li3() - (1.0_f64/z).li3() + (-z).cln().powi(3)/6.0_f64 + pi*pi/6.0_f64*(-z).cln()
    }
}


fn id3(z: Complex<f64>) -> Complex<f64> {
    if (1.0_f64 - z).re.abs() < std::f64::EPSILON || (z.re <= 0.0_f64 && z.im == 0.0_f64) {
        Complex::new(0.0_f64, 0.0_f64)
    } else {
        let pi = std::f64::consts::PI;
        let z3 = 1.202056903159594_f64;

        z.li3() + (1.0_f64 - z).li3() + (1.0_f64 - 1.0_f64/z).li3()
            - (z3 + z.cln().powi(3)/6.0_f64 + pi*pi/6.0_f64*z.cln()
               - 0.5_f64*z.cln().powi(2)*(1.0_f64 - z).cln())
    }
}


#[test]
fn special_values() {
    use num::Zero;
    let pi  = std::f64::consts::PI;
    let pi2 = pi*pi;
    let eps = 1e-15_f64;
    let ln2 = 2.0_f64.ln();
    let z3  = 1.202056903159594_f64;
    let phi = 0.5_f64*(5.0_f64.sqrt() + 1.0_f64); // golden ratio
    let zero = Complex::zero();

    assert_eq_complex!(zero.li3(), zero, eps);
    assert_eq_complex!(Complex::new(1.0_f64, 0.0_f64).li3(),
                       Complex::new(z3, 0.0_f64), eps);
    assert_eq_complex!(Complex::new(-1.0_f64, 0.0_f64).li3(),
                       Complex::new(-3.0_f64/4.0_f64*z3, 0.0_f64), eps);
    assert_eq_complex!(Complex::new(0.5_f64, 0.0_f64).li3(),
                       Complex::new(ln2.powi(3)/6.0_f64 - pi2/12.0_f64*ln2 + 7.0_f64/8.0_f64*z3, 0.0_f64), eps);
    assert_eq_complex!(Complex::new(1.0_f64/(phi*phi), 0.0_f64).li3(),
                       Complex::new(4.0_f64/5.0_f64*z3 + 2.0_f64/3.0_f64*phi.ln().powi(3) - 2.0_f64/15.0_f64*pi2*phi.ln(), 0.0_f64), eps);

    // test value that causes overflow if squared
    assert!(!Complex::new(1e300_f64, 1.0_f64).li3().is_infinite());
    assert!(!Complex::new(1.0_f64, 1e300_f64).li3().is_infinite());
    assert_eq_complex!(Complex::new(1e300_f64, 1.0_f64).li3(), Complex::new(-5.4934049431527088e7_f64, 749538.186928224_f64), eps);
    assert_eq_complex!(Complex::new(1.0_f64, 1e300_f64).li3(), Complex::new(-5.4936606061973454e7_f64, 374771.031356405_f64), eps);
}


#[test]
fn test_values() {
    let eps = 1e-14_f64;
    let values = common::read_data_file("Li3.txt").unwrap();

    for &(v, li3) in values.iter() {
        assert_eq_complex!(v.li3(), li3, eps);

        if v.im == 0.0_f64 {
            assert_eq_float!(v.re.li3(), li3.re, eps);
        }
    }
}


#[test]
fn identities() {
    use num::Zero;
    let eps = 1e-9_f64;
    let zero = Complex::<f64>::zero();
    let values = common::read_data_file("Li3.txt").unwrap();

    for &(v1, v2) in &values {
        assert_eq_complex!(id1(v1), zero, eps);
        assert_eq_complex!(id1(v2), zero, eps);
        assert_eq_complex!(id2(v1), zero, eps);
        assert_eq_complex!(id2(v2), zero, eps);
        assert_eq_complex!(id3(v1), zero, eps);
        assert_eq_complex!(id3(v2), zero, eps);
    }
}


#[test]
fn test_signed_zero() {
    let pz64 = 0.0_f64;
    let nz64 = -0.0_f64;

    assert!(pz64.li3().is_sign_positive());
    assert!(nz64.li3().is_sign_negative());

    assert!(Complex::new(pz64, pz64).li3().re.is_sign_positive());
    assert!(Complex::new(pz64, pz64).li3().im.is_sign_positive());
    assert!(Complex::new(pz64, nz64).li3().re.is_sign_positive());
    assert!(Complex::new(pz64, nz64).li3().im.is_sign_negative());
    assert!(Complex::new(nz64, pz64).li3().re.is_sign_negative());
    assert!(Complex::new(nz64, pz64).li3().im.is_sign_positive());
    assert!(Complex::new(nz64, nz64).li3().re.is_sign_negative());
    assert!(Complex::new(nz64, nz64).li3().im.is_sign_negative());
}
