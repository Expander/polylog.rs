extern crate num;
extern crate polylog;
extern crate rand;

use num::complex::Complex;
use polylog::{Li2, Li3, Li4};
use rand::Rng;
use std::time::{Duration, Instant};


#[test]
fn bench_real_li2() {
    let n = 1000000;
    let numbers = gen_real_numbers(-10.0, 10.0, n);

    let time: f64 = time(|| { let _: Vec<f64> = numbers.iter().map(|z| z.li2()).collect(); });

    println!("Evaluation of real Li2 {} times took: {}s", n, time);
}


#[test]
fn bench_complex_li2() {
    let n = 1000000;
    let numbers = gen_complex_numbers(-10.0, 10.0, n);

    let time: f64 = time(|| { let _: Vec<Complex<f64>> = numbers.iter().map(|z| z.li2()).collect(); });

    println!("Evaluation of complex Li2 {} times took: {}s", n, time);
}


#[test]
fn bench_complex_li3() {
    let n = 1000000;
    let numbers = gen_complex_numbers(-10.0, 10.0, n);

    let time: f64 = time(|| { let _: Vec<Complex<f64>> = numbers.iter().map(|z| z.li3()).collect(); });

    println!("Evaluation of complex Li3 {} times took: {}s", n, time);
}


#[test]
fn bench_complex_li4() {
    let n = 1000000;
    let numbers = gen_complex_numbers(-10.0, 10.0, n);

    let time: f64 = time(|| { let _: Vec<Complex<f64>> = numbers.iter().map(|z| z.li4()).collect(); });

    println!("Evaluation of complex Li4 {} times took: {}s", n, time);
}


fn gen_real_numbers(min: f64, max: f64, n: u32) -> Vec<f64> {
    let mut rng = rand::thread_rng();

    (0..n).map(|_| rng.gen_range(min, max)).collect()
}


fn gen_complex_numbers(min: f64, max: f64, n: u32) -> Vec<Complex<f64>> {
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
