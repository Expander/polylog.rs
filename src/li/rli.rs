use num::complex::Complex;
use crate::cln::CLn;
use crate::{Li0, Li1, Li2, Li3, Li4};
use super::eta::neg_eta;
use super::fac::{fac, inv_fac};
use super::harmonic::harmonic;
use super::zeta::zeta;

/// returns real n-th order polylogarithm Re[Li(n,x)] for real x
pub fn rli(n: i32, x: f64) -> f64 {
    let odd_sgn = |n| if is_even(n) { -1.0 } else { 1.0 };

    if x == 0.0 {
        0.0
    } else if x == 1.0 {
        zeta(n)
    } else if x == -1.0 {
        neg_eta(n)
    } else if x.is_nan() {
        std::f64::NAN
    } else if n < -1 {
        // arXiv:2010.09860
        let c = 4.0*std::f64::consts::PI*std::f64::consts::PI;
        let l2 = ln_sqr(x);
        if c*x*x < l2 {
            li_series(n, x)
        } else if l2 < 0.512*0.512*c {
            li_unity_neg(n, Complex::new(x, 0.0)).re
        } else {
            odd_sgn(n)*li_series(n, x.recip())
        }
    } else if n == -1 {
        x/((1.0 - x)*(1.0 - x))
    } else if n == 0 {
        x.li0()
    } else if n == 1 {
        x.li1()
    } else if n == 2 {
        x.li2()
    } else if n == 3 {
        x.li3()
    } else if n == 4 {
        x.li4()
    } else {
        // transform x to y in [-1,1]
        let (y, rest, sgn) = if x < -1.0 {
            (x.recip(), li_neg_rest(n, x), odd_sgn(n))
        } else if x < 1.0 {
            (x, 0.0, 1.0)
        } else { // x > 1.0
            (x.recip(), li_pos_rest(n, x), odd_sgn(n))
        };

        let li = if n < 20 && y > 0.75 {
            li_unity_pos(n, y)
        } else {
            li_series(n, y)
        };

        rest + sgn*li
    }
}

/// returns true if x is even, false otherwise
fn is_even(x: i32) -> bool {
    x & 1 == 0
}

/// returns |ln(x)|^2 for all x
fn ln_sqr(x: f64) -> f64 {
    if x < 0.0 {
        let l = (-x).ln();
        l*l + std::f64::consts::PI*std::f64::consts::PI
    } else if x == 0.0 {
        std::f64::NAN
    } else {
        let l = x.ln();
        l*l
    }
}

/// returns r.h.s. of inversion formula for x < -1:
///
/// Li(n,-x) + (-1)^n Li(n,-1/x)
///    = -ln(n,x)^n/n! + 2 sum(r=1:(n÷2), ln(x)^(n-2r)/(n-2r)! Li(2r,-1))
fn li_neg_rest(n: i32, x: f64) -> f64 {
    let l = (-x).ln();
    let l2 = l*l;

    if is_even(n) {
        let mut sum = 0.0;
        let mut p = 1.0; // collects l^(2u)
        for u in 0..n/2 {
            let old_sum = sum;
            sum += p*inv_fac(2*u)*neg_eta(n - 2*u);
            if sum == old_sum { break; }
            p *= l2;
        }
        2.0*sum - p*inv_fac(n)
    } else {
        let mut sum = 0.0;
        let mut p = l; // collects l^(2u + 1)
        for u in 0..(n - 1)/2 {
            let old_sum = sum;
            sum += p*inv_fac(2*u + 1)*neg_eta(n - 1 - 2*u);
            if sum == old_sum { break; }
            p *= l2;
        }
        2.0*sum - p*inv_fac(n)
    }
}

/// returns (sin((2n+1)x), cos((2n+1)x)), given
/// (sn, cn) = (sin(2nx), cos(2nx))   (previous value)
/// (s2, c2) = (sin(2x), sin(2x))     (initial value)
fn next_cosi((sn, cn): (f64, f64), (s2, c2): (f64, f64)) -> (f64, f64) {
    (sn*c2 + cn*s2, cn*c2 - sn*s2)
}

/// returns r.h.s. of inversion formula for x > 1;
/// same expression as in li_neg_rest(n,x), but with
/// complex logarithm ln(-x)
fn li_pos_rest(n: i32, x: f64) -> f64 {
    let pi = std::f64::consts::PI;
    let l = x.ln();
    let mag = l.hypot(pi); // |ln(-x)|
    let arg = pi.atan2(l); // arg(ln(-x))
    let l2 = mag*mag;      // |ln(-x)|^2

    if is_even(n) {
        let mut sum = 0.0;
        let mut p = 1.0; // collects mag^(2u)
        let mut cosi = (0.0, 1.0); // collects (sin(2*u*arg), cos(2*u*arg))
        let cosi2 = (2.0*arg).sin_cos();
        for u in 0..n/2 {
            let old_sum = sum;
            sum += p*cosi.1*inv_fac(2*u)*neg_eta(n - 2*u);
            if sum == old_sum { break; }
            p *= l2;
            cosi = next_cosi(cosi, cosi2);
        }
        2.0*sum - p*cosi.1*inv_fac(n)
    } else {
        let mut sum = 0.0;
        let mut p = mag; // collects mag^(2u + 1)
        let (s, c) = arg.sin_cos();
        let mut cosi = (s, c); // collects (sin((2*u + 1)*arg), cos((2*u + 1)*arg))
        let cosi2 = (2.0*s*c, 2.0*c*c - 1.0); // (2.0*arg).sin_cos()
        for u in 0..(n - 1)/2 {
            let old_sum = sum;
            sum += p*cosi.1*inv_fac(2*u + 1)*neg_eta(n - 1 - 2*u);
            if sum == old_sum { break; }
            p *= l2;
            cosi = next_cosi(cosi, cosi2);
        }
        2.0*sum - p*cosi.1*inv_fac(n)
    }
}

/// returns Li(n,x) using the series expansion of Li(n,x) for n > 0
/// and real x ~ 1 where 0 < x < 1:
///
/// Li(n,x) = sum(j=0:Inf, zeta(n-j) ln(x)^j/j!)
///
/// where
///
/// zeta(1) = -ln(-ln(x)) + harmonic(n - 1)
///
/// harmonic(n) = sum(k=1:n, 1/k)
fn li_unity_pos(n: i32, x: f64) -> f64 {
    let l = x.ln();
    let mut sum = zeta(n);
    let mut p = 1.0; // collects l^j/j!

    for j in 1..(n - 1) {
        p *= l/(j as f64);
        sum += zeta(n - j)*p;
    }

    p *= l/((n - 1) as f64);
    sum += (harmonic(n - 1) - (-l).ln())*p;

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
fn li_series(n: i32, x: f64) -> f64
{
    let mut sum = x;
    let mut xn = x*x;

    for k in 2..i32::MAX {
        let term = xn/(k as f64).powi(n);
        if !term.is_finite() { break; }
        let old_sum = sum;
        sum += term;
        if sum == old_sum { break; }
        xn *= x;
    }

    sum
}
