extern crate polylog;
extern crate num;
use polylog::Li;
mod common;

#[test]
fn test_values() {
    let eps = 1e-14;

    for n in 0..=4 {
        let filename = format!("Li{}.txt", n);
        let values = common::read_data_file(&filename).unwrap();

        for &(v, res) in values.iter() {
            if v.im == 0.0 {
                assert_eq_float!(v.re.li(n), res.re, eps);
            }
        }
    }
}
