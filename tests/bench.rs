extern crate num;
extern crate polylog;
extern crate rand;

use num::complex::Complex;
use polylog::{Li2, Li3, Li4, Li5, Li6};
use rand::Rng;
use std::time::{Duration, Instant};


#[test]
fn bench_real_li2() {
    let sample = gen_real_numbers(-10.0, 10.0, 1000000);
    bench_fn(|z: &f64| z.li2(), String::from("real Li2"), sample);
}


#[test]
fn bench_complex_li2() {
    let sample = gen_complex_numbers(-10.0, 10.0, 1000000);
    bench_fn(|z: &Complex<f64>| z.li2(), String::from("complex Li2"), sample);
}


#[test]
fn bench_complex_li3() {
    let sample = gen_complex_numbers(-10.0, 10.0, 1000000);
    bench_fn(|z: &Complex<f64>| z.li3(), String::from("complex Li3"), sample);
}


#[test]
fn bench_complex_li4() {
    let sample = gen_complex_numbers(-10.0, 10.0, 1000000);
    bench_fn(|z: &Complex<f64>| z.li4(), String::from("complex Li4"), sample);
}


#[test]
fn bench_complex_li5() {
    let sample = gen_complex_numbers(-10.0, 10.0, 1000000);
    bench_fn(|z: &Complex<f64>| z.li5(), String::from("complex Li5"), sample);
}


#[test]
fn bench_complex_li6() {
    let sample = gen_complex_numbers(-10.0, 10.0, 1000000);
    bench_fn(|z: &Complex<f64>| z.li6(), String::from("complex Li6"), sample);
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


fn bench_fn<F, S>(fun: F, name: String, sample: Vec<S>)
    where F: Fn(&S) -> S
{
    let sample_size = sample.len();
    let time: f64 = time(|| { let _: Vec<S> = sample.iter().map(|z| fun(z)).collect(); });
    println!("Evaluation of {} {} times took: {}s", name, sample_size, time); 
}
