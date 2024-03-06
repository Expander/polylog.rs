use num::complex::Complex;
use polylog::Li;
mod common;

#[test]
fn test_values() {
    struct Ni { n: i32, eps: f64 }

    let ni = vec![
        Ni { n:     -10, eps: 1e-09_f64},
        Ni { n:      -9, eps: 1e-10_f64},
        Ni { n:      -8, eps: 1e-10_f64},
        Ni { n:      -7, eps: 1e-12_f64},
        Ni { n:      -6, eps: 1e-12_f64},
        Ni { n:      -5, eps: 1e-10_f64},
        Ni { n:      -4, eps: 1e-13_f64},
        Ni { n:      -3, eps: 1e-13_f64},
        Ni { n:      -2, eps: 1e-13_f64},
        Ni { n:      -1, eps: 1e-14_f64},
        Ni { n:       0, eps: 1e-14_f64},
        Ni { n:       1, eps: 1e-14_f64},
        Ni { n:       2, eps: 1e-14_f64},
        Ni { n:       3, eps: 1e-14_f64},
        Ni { n:       4, eps: 1e-14_f64},
        Ni { n:       5, eps: 1e-14_f64},
        Ni { n:       6, eps: 1e-14_f64},
        Ni { n:     100, eps: 1e-14_f64},
        Ni { n: 1000000, eps: 1e-14_f64},
    ];

    for n in ni.into_iter() {
        let filename = format!("Li{}.txt", n.n);
        let values = common::read_data_file(&filename).unwrap();

        for &(v, res) in values.iter() {
            assert_eq_complex!(v.li(n.n), res, n.eps);

            if v.im == 0.0_f64 {
                assert_eq_float!(v.re.li(n.n), res.re, n.eps);
            }
        }

        assert!(std::f64::NAN.li(n.n).is_nan());
    }

    // value close to boundary between series 1 and 2 in arXiv:2010.09860
    assert_eq_float!((-0.50001_f64).li(-2), -0.074072592582716422_f64, 1e-14_f64);
    assert_eq_complex!(Complex::<f64>::new(-0.50001_f64, 0.0_f64).li(-2),
                       Complex::<f64>::new(-0.074072592582716422_f64, 0.0_f64), 1e-14_f64);

    // value sensitive to proper treatment of 0.0 vs -0.0 in imag(z)
    let z = Complex::new(1.5_f64, 0.0_f64);
    assert_eq_complex!(z.li(10), Complex::<f64>::new(1.5022603281703005298_f64, -2.56429642116111388671e-9_f64), 1e-14_f64);
    assert_eq_complex!((-z).li(10), Complex::<f64>::new(-1.4978556954869267594_f64, 0.0_f64), 1e-14_f64);

    // test value that causes overflow if squared
    assert!(!Complex::new(1e300_f64, 1.0_f64).li(7).is_infinite());
    assert!(!Complex::new(1.0_f64, 1e300_f64).li(7).is_infinite());
    assert_eq_complex!(Complex::new(1e300_f64, 1.0_f64).li(7), Complex::new(-1.4886831990993457e16_f64, 4.74066248802866e14_f64), 1e-15_f64);
    assert_eq_complex!(Complex::new(1.0_f64, 1e300_f64).li(7), Complex::new(-1.489168315226607e16_f64, 2.3705150998401e14_f64), 1e-5_f64);

    // test non-finite input
    assert!(Complex::new(f64::NAN, f64::NAN).li(7).is_nan());
    assert!(Complex::new(f64::INFINITY, f64::INFINITY).li(7).is_infinite());
}


#[test]
fn test_signed_zero() {
    let pz64 = 0.0_f64;
    let nz64 = -0.0_f64;

    for n in (-100..100).into_iter() {
        assert!(pz64.li(n).is_sign_positive());
        assert!(nz64.li(n).is_sign_negative());

        assert!(Complex::new(pz64, pz64).li(n).re.is_sign_positive());
        assert!(Complex::new(pz64, pz64).li(n).im.is_sign_positive());
        assert!(Complex::new(pz64, nz64).li(n).re.is_sign_positive());
        assert!(Complex::new(pz64, nz64).li(n).im.is_sign_negative());
        assert!(Complex::new(nz64, pz64).li(n).re.is_sign_negative());
        assert!(Complex::new(nz64, pz64).li(n).im.is_sign_positive());
        assert!(Complex::new(nz64, nz64).li(n).re.is_sign_negative());
        assert!(Complex::new(nz64, nz64).li(n).im.is_sign_negative());
    }
}
