use num::complex::Complex;
use polylog::{Li2, Li3, Li4, Li5, Li6, Li};
use criterion::*;


fn bench_real_li2(c: &mut Criterion) {
    let mut group = c.benchmark_group("li2(x)");
    group.bench_function("x=0.25", |b| b.iter(|| black_box(0.25).li2()));
    group.finish();
}


fn bench_complex_li2(c: &mut Criterion) {
    let mut group = c.benchmark_group("li2(z)");
    group.bench_function("z=0.25+0.25i", |b| b.iter(|| black_box(Complex::new(0.25, 0.25)).li2()));
    group.bench_function("z=-0.7+0.7i" , |b| b.iter(|| black_box(Complex::new(-0.7,  0.7)).li2()));
    group.finish();
}


fn bench_real_li3(c: &mut Criterion) {
    let mut group = c.benchmark_group("li3(x)");
    group.bench_function("x=0.25", |b| b.iter(|| black_box(0.25).li3()));
    group.bench_function("x=-0.5", |b| b.iter(|| black_box(-0.5).li3()));
    group.finish();
}


fn bench_complex_li3(c: &mut Criterion) {
    let mut group = c.benchmark_group("li3(z)");
    group.bench_function("z=0.25+0.25i", |b| b.iter(|| black_box(Complex::new(0.25, 0.25)).li3()));
    group.bench_function("z=-0.7+0.7i" , |b| b.iter(|| black_box(Complex::new(-0.7,  0.7)).li3()));
    group.finish();
}


fn bench_real_li4(c: &mut Criterion) {
    let mut group = c.benchmark_group("li4(x)");
    group.bench_function("x=0.5" , |b| b.iter(|| black_box(0.25).li4()));
    group.bench_function("x=-0.5", |b| b.iter(|| black_box(-0.5).li4()));
    group.finish();
}


fn bench_complex_li4(c: &mut Criterion) {
    let mut group = c.benchmark_group("li4(z)");
    group.bench_function("z=0.25+0.25i", |b| b.iter(|| black_box(Complex::new(0.25, 0.25)).li4()));
    group.bench_function("z=-0.7+0.7i" , |b| b.iter(|| black_box(Complex::new(-0.7,  0.7)).li4()));
    group.finish();
}


fn bench_complex_li5(c: &mut Criterion) {
    let mut group = c.benchmark_group("li5(z)");
    group.bench_function("z=0.25+0.25i", |b| b.iter(|| black_box(Complex::new(0.25, 0.25)).li5()));
    group.bench_function("z=-0.7+0.7i" , |b| b.iter(|| black_box(Complex::new(-0.7,  0.7)).li5()));
    group.finish();
}


fn bench_complex_li6(c: &mut Criterion) {
    let mut group = c.benchmark_group("li6(z)");
    group.bench_function("z=0.25+0.25i", |b| b.iter(|| black_box(Complex::new(0.25, 0.25)).li6()));
    group.bench_function("z=-0.7+0.7i" , |b| b.iter(|| black_box(Complex::new(-0.7,  0.7)).li6()));
    group.finish();
}


fn bench_real_li(c: &mut Criterion) {

    let mut ni: Vec<_> = (-10..10).step_by(2).collect();
    let n2 = vec![-10000, -1000, -100, 100, 1000, 1000_000];
    ni.extend(n2);

    let mut group = c.benchmark_group("li(n,x)");

    for n in ni.into_iter() {
        group.bench_function(format!("n={},x=0.2", n), |b| b.iter(|| black_box(0.2).li(n)));
        group.bench_function(format!("n={},x=0.8", n), |b| b.iter(|| black_box(0.8).li(n)));
    }

    group.finish();
}


fn bench_complex_li(c: &mut Criterion) {

    let mut ni: Vec<_> = (-10..10).step_by(2).collect();
    let n2 = vec![-10000, -1000, -100, 100, 1000, 1000_000];
    ni.extend(n2);

    let mut group = c.benchmark_group("li(n,z)");

    for n in ni.into_iter() {
        group.bench_function(format!("n={},z=0.25+0.25i", n), |b| b.iter(|| black_box(Complex::new(0.25, 0.25)).li(n)));
        group.bench_function(format!("n={},z=-0.7+0.7i" , n), |b| b.iter(|| black_box(Complex::new(-0.7,  0.7)).li(n)));
    }

    group.finish();
}


criterion_group!(benches,
                 bench_real_li2, bench_complex_li2,
                 bench_real_li3, bench_complex_li3,
                 bench_real_li4, bench_complex_li4,
                 bench_complex_li5,
                 bench_complex_li6,
                 bench_real_li,
                 bench_complex_li
);
criterion_main!(benches);
