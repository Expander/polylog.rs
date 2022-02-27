extern crate polylog;
extern crate num;
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
            if n.n <= 6 || v.norm_sqr() <= 0.75*0.75 {
                // println!("n = {}, z = {}", n.n, v);
                assert_eq_complex!(v.li(n.n), res, n.eps);
            }
        }

        assert!(std::f64::NAN.li(n.n).is_nan());
    }

    // value close to boundary between series 1 and 2 in arXiv:2010.09860
    assert_eq_float!((-0.50001_f64).li(-2), -0.074072592582716422_f64, 1e-14);
}
