use num::complex::Complex;
use crate::cln::CLn;
use crate::{Li0, Li1, Li2, Li3, Li4, Li5, Li6};
use super::eta::neg_eta;
use super::fac::{fac, inv_fac};
use super::harmonic::harmonic;
use super::zeta::zeta;

/// returns complex n-th order polylogarithm Li(n,z) for complex z
pub fn cli(n: i32, z: Complex<f64>) -> Complex<f64> {
    if z.is_nan() {
        Complex::new(f64::NAN, f64::NAN)
    } else if z.is_infinite() {
        Complex::new(f64::NEG_INFINITY, 0.0)
    } else if z == Complex::new(0.0, 0.0) {
        Complex::new(0.0, 0.0)
    } else if z == Complex::new(1.0, 0.0) {
        if n <= 0 {
            Complex::new(f64::INFINITY, f64::INFINITY)
        } else {
            Complex::new(zeta(n), 0.0)
        }
    } else if z == Complex::new(-1.0, 0.0) {
        Complex::new(neg_eta(n), 0.0)
    } else if n < -1 {
        // arXiv:2010.09860
        let c = 4.0*std::f64::consts::PI*std::f64::consts::PI;
        let l2 = z.cln().norm_sqr();
        if c*z.norm_sqr() < l2 {
            li_series(n, z)
        } else if l2 < 0.512*0.512*c {
            li_unity_neg(n, z)
        } else {
            let sqrtz = z.sqrt();
            2.0_f64.powi(n - 1)*(cli(n, sqrtz) + cli(n, -sqrtz))
        }
    } else if n == -1 {
        z/((1.0 - z)*(1.0 - z))
    } else if n == 0 {
        z.li0()
    } else if n == 1 {
        z.li1()
    } else if n == 2 {
        z.li2()
    } else if n == 3 {
        z.li3()
    } else if n == 4 {
        z.li4()
    } else if n == 5 {
        z.li5()
    } else if n == 6 {
        z.li6()
    } else if z.norm_sqr() <= 0.75*0.75 {
        li_series(n, z)
    } else if z.norm_sqr() >= 1.4*1.4 {
        let sgn = if is_even(n) { -1.0 } else { 1.0 };
        sgn*li_series(n, 1.0/z) + li_rest(n, z)
    } else {
        li_unity_pos(n, z)
    }
}

/// returns true if x is even, false otherwise
fn is_even(x: i32) -> bool {
    x & 1 == 0
}

/// returns r.h.s. of inversion formula for complex z
///
/// Li(n,-z) + (-1)^n Li(n,-1/z)
///    = -ln(n,z)^n/n! + 2 sum(k=1:(n÷2), ln(z)^(n-2k)/(n-2k)! Li(2k,-1))
fn li_rest(n: i32, z: Complex<f64>) -> Complex<f64> {
    let lnz = (-z).cln();
    let lnz2 = lnz*lnz;
    let kmax = if is_even(n) { n/2 } else { (n - 1)/2 };
    let mut p = if is_even(n) { Complex::new(1.0, 0.0) } else { lnz };
    let mut sum = Complex::new(0.0, 0.0);

    for k in (1..=kmax).rev() {
        let ifac = inv_fac(n - 2*k);
        if ifac == 0.0 { return 2.0*sum; }
        sum += neg_eta(2*k)*ifac*p;
        p *= lnz2;
    }

    2.0*sum - p*inv_fac(n)
}

/// returns Li(n,z) using the series expansion of Li(n,z) for n > 0
/// and complex z ~ 1:
///
/// Li(n,z) = sum(j=0:Inf, zeta(n-j) ln(z)^j/j!)
///
/// where
///
/// zeta(1) = -ln(-ln(z)) + harmonic(n - 1)
///
/// harmonic(n) = sum(k=1:n, 1/k)
fn li_unity_pos(n: i32, z: Complex<f64>) -> Complex<f64> {
    let l = z.cln();
    let mut sum = Complex::new(zeta(n), 0.0);
    let mut p = Complex::new(1.0, 0.0); // collects l^j/j!

    for j in 1..(n - 1) {
        p *= l/(j as f64);
        let old_sum = sum;
        sum += zeta(n - j)*p;
        if sum == old_sum { break; }
    }

    p *= l/((n - 1) as f64);
    sum += (harmonic(n - 1) - (-l).cln())*p;

    p *= l/(n as f64);
    sum += zeta(0)*p;

    p *= l/((n + 1) as f64);
    sum += zeta(-1)*p;

    let l2 = l*l;

    for j in ((n + 3)..i32::MAX).step_by(2) {
        p *= l2/(((j - 1)*j) as f64);
        let old_sum = sum;
        sum += zeta(n - j)*p;
        if sum == old_sum { break; }
    }

    sum
}

/// returns Li(n,x) using the series expansion for n < 0 and x ~ 1
///
/// Li(n,x) = gamma(1-n) (-ln(x))^(n-1)
///           + sum(k=0:Inf, zeta(n-k) ln(x)^k/k!)
fn li_unity_neg(n: i32, z: Complex<f64>) -> Complex<f64> {
    let lnz = z.cln();
    let lnz2 = lnz*lnz;
    let mut sum = fac(-n)*(-lnz).powi(n - 1);
    let (mut k, mut lnzk) = if is_even(n) {
        (1, lnz)
    } else {
        sum += zeta(n);
        (2, lnz2)
    };

    loop {
        let term = zeta(n - k)*inv_fac(k)*lnzk;
        if !term.is_finite() { break; }
        let sum_old = sum;
        sum += term;
        if sum == sum_old || k >= i32::MAX - 2 { break; }
        lnzk *= lnz2;
        k += 2;
    }

    sum
}

/// returns Li(n,x) using the naive series expansion of Li(n,x)
/// for |x| < 1:
///
/// Li(n,x) = sum(k=1:Inf, x^k/k^n)
fn li_series(n: i32, z: Complex<f64>) -> Complex<f64>
{
    let mut sum = z;
    let mut zn = z*z;

    for k in 2..i32::MAX {
        let term = zn/(k as f64).powi(n);
        if !term.is_finite() { break; }
        let old_sum = sum;
        sum += term;
        if sum == old_sum { break; }
        zn *= z;
    }

    sum
}
