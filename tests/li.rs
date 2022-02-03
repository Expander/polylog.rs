extern crate polylog;
extern crate num;
use polylog::Li;
mod common;

#[test]
fn test_values() {
    let eps = 1e-14;
    let ni = vec![0, 1, 2, 3, 4, 5, 6, 100];

    for n in ni.into_iter() {
        let filename = format!("Li{}.txt", n);
        let values = common::read_data_file(&filename).unwrap();

        for &(v, res) in values.iter() {
            if v.im == 0.0 {
                if n <= 4 || (v.re <= 1.0 && v.re >= -1.0) {
                    assert_eq_float!(v.re.li(n), res.re, eps);
                }
            }
        }

        assert!(std::f64::NAN.li(n).is_nan());
    }
}
