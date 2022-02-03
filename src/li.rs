use {Li0, Li1, Li2, Li3, Li4};
use {harmonic, zeta};

/// Provides the n-th order polylogarithm function `li()` of a number of type `T`.
pub trait Li<T> {
    fn li(&self, n: i32) -> T;
}

impl Li<f64> for f64 {
    /// Returns the real n-th order polylogarithm of a real number of
    /// type `f64`.
    ///
    /// # Example:
    /// ```
    /// use polylog::Li;
    ///
    /// let z = 1.0;
    /// let n = 10;
    /// println!("Li({},{}) = {}", n, z, z.li(n));
    /// ```
    fn li(&self, n: i32) -> f64 {
        if n < 0 {
            panic!("li(n) not implemented for n < 0 (given value: n = {})", n);
        } else if n == 0 {
            self.li0()
        } else if n == 1 {
            self.li1()
        } else if n == 2 {
            self.li2()
        } else if n == 3 {
            self.li3()
        } else if n == 4 {
            self.li4()
        } else if *self == 1.0 {
            zeta::zeta(n)
        } else if *self == -1.0 {
            li_minus_1(n)
        } else if self.is_nan() {
            std::f64::NAN
        } else {
            let x = *self;
            let is_odd = |x: i32| x & 1 == 1;

            // transform x to y in [-1,1]
            let (y, rest, sgn) = if x < -1.0 {
                (1.0/x, li_neg_rest(n, x), if is_odd(n) { 1.0 } else { -1.0 })
            } else if x < 1.0 {
                (x, 0.0, 1.0)
            } else { // x > 1.0
                (1.0/x, li_pos_rest(n, x), if is_odd(n) { 1.0 } else { -1.0})
            };

            let li = if n < 20 && y > 0.75 {
                li_series_one(n, y)
            } else {
                li_series_naive(n, y)
            };

            rest + sgn*li
        }
    }
}

/// returns r.h.s. of inversion formula for x < -1:
///
/// Li(n,-x) + (-1)^n Li(n,-1/x)
///    = -log(n,x)^n/n! + 2 sum(r=1:(n÷2), log(x)^(n-2r)/(n-2r)! Li(2r,-1))
fn li_neg_rest(n: i32, x: f64) -> f64 {
    0.0
}

/// returns r.h.s. of inversion formula for x > 1;
/// same expression as in li_neg_rest(n,x), but with
/// complex logarithm log(Complex(-x))
fn li_pos_rest(n: i32, x: f64) -> f64 {
    0.0
}

/// returns Li(n,x) using the series expansion of Li(n,x) for x ~ 1
/// where 0 < x < 1:
///
/// Li(n,x) = sum(j=0:Inf, zeta(n-j) log(x)^j/j!)
///
/// where
///
/// zeta(1) = -log(-log(x)) + harmonic(n - 1)
///
/// harmonic(n) = sum(k=1:n, 1/k)
fn li_series_one(n: i32, x: f64) -> f64 {
    let l = x.ln();
    let mut sum = zeta::zeta(n);
    let mut p = 1.0; // collects l^j/j!

    for j in 1..(n - 1) {
        p *= l/(j as f64);
        sum += zeta::zeta(n - j)*p;
    }

    p *= l/((n - 1) as f64);
    sum += (harmonic::harmonic(n - 1) - (-l).ln())*p;

    p *= l/(n as f64);
    sum += zeta::zeta(0)*p;

    p *= l/((n + 1) as f64);
    sum += zeta::zeta(-1)*p;

    let l2 = l*l;
    let mut old_sum;

    for j in ((n + 3)..i32::MAX).step_by(2) {
        p *= l2/(((j - 1)*j) as f64);
        old_sum = sum;
        sum += zeta::zeta(n - j)*p;
        if sum == old_sum {
            break;
        }
    }

    sum
}

/// returns Li(n,x) using the naive series expansion of Li(n,x)
/// for |x| < 1:
///
/// Li(n,x) = sum(k=1:Inf, x^k/k^n)
fn li_series_naive(n: i32, x: f64) -> f64 {
    let mut sum = x;
    let mut old_sum;
    let mut xn = x*x;

    for k in 2..i32::MAX {
        old_sum = sum;
        sum += xn/(k as f64).powi(n);
        if sum == old_sum {
            break;
        }
        xn *= x;
    }

    sum
}

// Table[PolyLog[n,-1], {n,1,54}]
const LI_MINUS_1_COEFF: [f64; 54] = [
    -0.69314718055994531, -0.82246703342411322, -0.90154267736969571,
    -0.94703282949724592, -0.97211977044690931, -0.98555109129743510,
    -0.99259381992283028, -0.99623300185264790, -0.99809429754160533,
    -0.99903950759827157, -0.99951714349806075, -0.99975768514385819,
    -0.99987854276326512, -0.99993917034597972, -0.99996955121309924,
    -0.99998476421490611, -0.99999237829204101, -0.99999618786961011,
    -0.99999809350817168, -0.99999904661158152, -0.99999952325821554,
    -0.99999976161323082, -0.99999988080131844, -0.99999994039889239,
    -0.99999997019885696, -0.99999998509923200, -0.99999999254955048,
    -0.99999999627475340, -0.99999999813736942, -0.99999999906868228,
    -0.9999999995343403 , -0.9999999997671699 , -0.9999999998835849 ,
    -0.9999999999417924 , -0.9999999999708962 , -0.9999999999854481 ,
    -0.9999999999927240 , -0.9999999999963620 , -0.9999999999981810 ,
    -0.9999999999990905 , -0.9999999999995453 , -0.9999999999997726 ,
    -0.9999999999998863 , -0.9999999999999432 , -0.9999999999999716 ,
    -0.9999999999999858 , -0.9999999999999929 , -0.9999999999999964 ,
    -0.9999999999999982 , -0.9999999999999991 , -0.9999999999999996 ,
    -0.9999999999999998 , -0.9999999999999999 , -0.9999999999999999
];

// returns Li(n,-1) = (2.0^(1 - n) - 1.0)*zeta(n) for n > 0
fn li_minus_1(n: i32) -> f64 {
    if n < 1 {
        panic!("li_minus_1 not implemented for n < 1 (given value: n = {})", n);
    } else if n as usize <= LI_MINUS_1_COEFF.len() {
        LI_MINUS_1_COEFF[(n - 1) as usize]
    } else {
        -1.0
    }
}

#[test]
fn test_li_minus_1() {
    assert!(li_minus_1(1)  == -0.69314718055994531);
    assert!(li_minus_1(2)  == -0.82246703342411322);
    assert!(li_minus_1(52) == -0.9999999999999998);
    assert!(li_minus_1(53) == -0.9999999999999999);
    assert!(li_minus_1(54) == -0.9999999999999999);
    assert!(li_minus_1(55) == -1.0);
    assert!(li_minus_1(56) == -1.0);
}
