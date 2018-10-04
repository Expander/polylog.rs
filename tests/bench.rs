extern crate num;
extern crate polylog;
extern crate rand;

use num::complex::Complex;
use polylog::{Li2, Li3};
use rand::Rng;
use std::time::{Duration, Instant};


#[test]
fn bench_li2() {
    let n = 1000000;
    let numbers = gen_numbers(-10.0, 10.0, n);

    let time: f64 = numbers.iter().map(|z| time(|| { z.li2(); })).sum();

    println!("Evaluation of Li2 {} times took: {}ms", n, time);
}


#[test]
fn bench_li3() {
    let n = 1000000;
    let numbers = gen_numbers(-10.0, 10.0, n);

    let time: f64 = numbers.iter().map(|z| time(|| { z.li3(); })).sum();

    println!("Evaluation of Li3 {} times took: {}ms", n, time);
}


fn gen_numbers(min: f64, max: f64, n: u32) -> Vec<Complex<f64>> {
    let mut rng = rand::thread_rng();

    (0..n).map(|_| {
        Complex::new(
            rng.gen_range(min, max),
            rng.gen_range(min, max)
        )        
    }).collect()
}


fn duration_nanos(d: Duration) -> f64 {
    (d.as_secs() as f64) + f64::from(d.subsec_nanos()) / 1e9
}


fn time<F>(mut f: F) -> f64
    where F: FnMut()
{
    let start = Instant::now();
    f();
    let end = Instant::now();
    let dur = duration_nanos(end.duration_since(start));
    f64::from(dur)
}
