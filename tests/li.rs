extern crate polylog;
extern crate num;
use num::complex::Complex;
use num::Float;
use polylog::Li;
mod common;

#[test]
fn test_values() {
    struct Ni { n: i32, eps: f64 }

    let ni = vec![
        Ni { n: -10, eps: 1e-09},
        Ni { n:  -9, eps: 1e-10},
        Ni { n:  -8, eps: 1e-10},
        Ni { n:  -7, eps: 1e-12},
        Ni { n:  -6, eps: 1e-12},
        Ni { n:  -5, eps: 1e-10},
        Ni { n:  -4, eps: 1e-13},
        Ni { n:  -3, eps: 1e-13},
        Ni { n:  -2, eps: 1e-13},
        Ni { n:  -1, eps: 1e-14},
        Ni { n:   0, eps: 1e-14},
        Ni { n:   1, eps: 1e-14},
        Ni { n:   2, eps: 1e-14},
        Ni { n:   3, eps: 1e-14},
        Ni { n:   4, eps: 1e-14},
        Ni { n:   5, eps: 1e-14},
        Ni { n:   6, eps: 1e-14},
        Ni { n: 100, eps: 1e-14},
    ];

    for n in ni.into_iter() {
        let filename = format!("Li{}.txt", n.n);
        let values = common::read_data_file(&filename).unwrap();

        for &(v, res) in values.iter() {
            if v.im == 0.0 {
                assert_eq_float!(v.re.li(n.n), res.re, n.eps);
            }
            assert_eq_complex!(v.li(n.n), res, n.eps);
        }

        assert!(std::f64::NAN.li(n.n).is_nan());
    }

    // value close to boundary between series 1 and 2 in arXiv:2010.09860
    assert_eq_float!((-0.50001_f64).li(-2), -0.074072592582716422_f64, 1e-14);
    assert_eq_complex!(Complex::<f64>::new(-0.50001, 0.0).li(-2),
                       Complex::<f64>::new(-0.074072592582716422, 0.0), 1e-14);

    // value sensitive to proper treatment of 0.0 vs -0.0 in imag(z)
    let z = Complex::new(1.5, 0.0);
    assert_eq_complex!(z.li(10), Complex::<f64>::new(1.5022603281703005298, -2.56429642116111388671e-9), 1e-14);
    assert_eq_complex!((-z).li(10), Complex::<f64>::new(-1.4978556954869267594, 0.0), 1e-14);

    // test value that causes overflow if squared
    assert!(!Complex::new(1e300, 1.0).li(7).is_infinite());
    assert!(!Complex::new(1.0, 1e300).li(7).is_infinite());
    assert_eq_float!(Complex::new(1e300, 1.0).li(7).re, Complex::new(-1.4886831990993457e16, -4.74066248802866e14).re, 1e-15);
    assert_eq_float!(Complex::new(1.0, 1e300).li(7).re, Complex::new(-1.489168315226607e16, 2.3705150998401e14).re, 1e-5);
}
