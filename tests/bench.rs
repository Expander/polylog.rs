extern crate num;
extern crate polylog;
extern crate rand;

use num::complex::Complex;
use polylog::{Li2, Li3, Li4, Li5, Li6, Li};
use rand::Rng;
use std::time::{Duration, Instant};


#[test]
fn bench_real_li2() {
    let sample = gen_real_numbers(0.0, 0.5, 10_000_000);
    bench_fn(|z: &f64| z.li2(), String::from("real Li2"), sample);
}


#[test]
fn bench_complex_li2() {
    let sample = gen_complex_numbers(-1.0, 1.0, 10_000_000);
    bench_fn(|z: &Complex<f64>| z.li2(), String::from("complex Li2"), sample);
}


#[test]
fn bench_real_li3() {
    let sample = gen_real_numbers(-1.0, 0.5, 10_000_000);
    bench_fn(|z: &f64| z.li3(), String::from("real Li3"), sample);
}


#[test]
fn bench_complex_li3() {
    let sample = gen_complex_numbers(-1.0, 1.0, 10_000_000);
    bench_fn(|z: &Complex<f64>| z.li3(), String::from("complex Li3"), sample);
}


#[test]
fn bench_real_li4() {
    let sample = gen_real_numbers(-1.0, 1.0, 10_000_000);
    bench_fn(|z: &f64| z.li4(), String::from("real Li4"), sample);
}


#[test]
fn bench_complex_li4() {
    let sample = gen_complex_numbers(-1.0, 1.0, 10_000_000);
    bench_fn(|z: &Complex<f64>| z.li4(), String::from("complex Li4"), sample);
}


#[test]
fn bench_complex_li5() {
    let sample = gen_complex_numbers(-1.0, 1.0, 10_000_000);
    bench_fn(|z: &Complex<f64>| z.li5(), String::from("complex Li5"), sample);
}


#[test]
fn bench_complex_li6() {
    let sample = gen_complex_numbers(-1.0, 1.0, 10_000_000);
    bench_fn(|z: &Complex<f64>| z.li6(), String::from("complex Li6"), sample);
}


#[test]
fn bench_real_li() {
    println!("Benchmark of Li(n,x):");

    let mut ni: Vec<_> = (5..10).collect();
    let n2 = vec![100, 1000, 1000_000];
    ni.extend(n2);

    for n in ni.into_iter() {
        let sample = gen_real_numbers(-1.0, 1.0, 10_000_000);
        bench_fn(|z: &f64| z.li(n), format!("real Li_{}", n), sample);
    }
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
