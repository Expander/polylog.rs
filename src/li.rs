use std::ops::{AddAssign, Div, Mul, MulAssign};
use std::cmp::PartialEq;
use num::complex::Complex;
use crate::cln::CLn;
use crate::{Li0, Li1, Li2, Li3, Li4, Li5, Li6};
mod harmonic;
mod fac;
mod zeta;

/// Provides the n-th order polylogarithm function `li()` of a number of type `T`.
pub trait Li<T> {
    fn li(&self, n: i32) -> T;
}

impl Li<Complex<f64>> for Complex<f64> {
    /// Returns the complex n-th order polylogarithm of a complex
    /// number of type `Complex<f64>` for all integers `n`.
    ///
    /// The implementation is an adaption of [[arxiv:2010.09860]].
    ///
    /// [arxiv:2010.09860]: https://arxiv.org/abs/2010.09860
    ///
    /// # Example:
    /// ```
    /// extern crate num;
    /// use num::complex::Complex;
    /// use polylog::Li;
    ///
    /// let z = Complex::new(1.0, 1.0);
    /// let n = 10;
    /// println!("Li({},{}) = {}", n, z, z.li(n));
    /// ```
    fn li(&self, n: i32) -> Complex<f64> {
        if self.is_nan() {
            Complex::new(f64::NAN, f64::NAN)
        } else if self.is_infinite() {
            Complex::new(f64::NEG_INFINITY, 0.0)
        } else if *self == Complex::new(0.0, 0.0) {
            Complex::new(0.0, 0.0)
        } else if *self == Complex::new(1.0, 0.0) {
            if n <= 0 {
                Complex::new(f64::INFINITY, f64::INFINITY)
            } else {
                Complex::new(zeta::zeta(n), 0.0)
            }
        } else if *self == Complex::new(-1.0, 0.0) {
            Complex::new(li_minus_1(n), 0.0)
        } else if n < -1 {
            // arXiv:2010.09860
            let c = 4.0*std::f64::consts::PI*std::f64::consts::PI;
            let z = *self;
            let l2 = z.cln().norm_sqr();
            if c*z.norm_sqr() < l2 {
                li_series(n, z)
            } else if l2 < 0.512*0.512*c {
                li_unity_neg(n, z)
            } else {
                let sqrtz = z.sqrt();
                2.0_f64.powi(n - 1)*(sqrtz.li(n) + (-sqrtz).li(n))
            }
        } else if n == -1 {
            let z = *self;
            z/((1.0 - z)*(1.0 - z))
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
        } else if n == 5 {
            self.li5()
        } else if n == 6 {
            self.li6()
        } else if self.norm_sqr() <= 0.75*0.75 {
            li_series(n, *self)
        } else {
            Complex::new(0.0, 0.0) // @todo
        }
    }
}

impl Li<f64> for f64 {
    /// Returns the real n-th order polylogarithm of a real number of
    /// type `f64` for all integers `n`.
    ///
    /// The implementation for `n < 0` is an adaption of
    /// [[arxiv:2010.09860]].
    ///
    /// [arxiv:2010.09860]: https://arxiv.org/abs/2010.09860
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
        let odd_sgn = |n| if is_even(n) { -1.0 } else { 1.0 };

        if *self == 0.0 {
            0.0
        } else if *self == 1.0 {
            zeta::zeta(n)
        } else if *self == -1.0 {
            li_minus_1(n)
        } else if self.is_nan() {
            std::f64::NAN
        } else if n < -1 {
            // arXiv:2010.09860
            let x = *self;
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
            *self/((1.0 - *self)*(1.0 - *self))
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
        } else {
            let x = *self;

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
            sum += p*fac::inv_fac(2*u)*li_minus_1(n - 2*u);
            if sum == old_sum { break; }
            p *= l2;
        }
        2.0*sum - p*fac::inv_fac(n)
    } else {
        let mut sum = 0.0;
        let mut p = l; // collects l^(2u + 1)
        for u in 0..(n - 1)/2 {
            let old_sum = sum;
            sum += p*fac::inv_fac(2*u + 1)*li_minus_1(n - 1 - 2*u);
            if sum == old_sum { break; }
            p *= l2;
        }
        2.0*sum - p*fac::inv_fac(n)
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
            sum += p*cosi.1*fac::inv_fac(2*u)*li_minus_1(n - 2*u);
            if sum == old_sum { break; }
            p *= l2;
            cosi = next_cosi(cosi, cosi2);
        }
        2.0*sum - p*cosi.1*fac::inv_fac(n)
    } else {
        let mut sum = 0.0;
        let mut p = mag; // collects mag^(2u + 1)
        let (s, c) = arg.sin_cos();
        let mut cosi = (s, c); // collects (sin((2*u + 1)*arg), cos((2*u + 1)*arg))
        let cosi2 = (2.0*s*c, 2.0*c*c - 1.0); // (2.0*arg).sin_cos()
        for u in 0..(n - 1)/2 {
            let old_sum = sum;
            sum += p*cosi.1*fac::inv_fac(2*u + 1)*li_minus_1(n - 1 - 2*u);
            if sum == old_sum { break; }
            p *= l2;
            cosi = next_cosi(cosi, cosi2);
        }
        2.0*sum - p*cosi.1*fac::inv_fac(n)
    }
}

/// returns Li(n,x) using the series expansion of Li(n,x) for n > 0
/// and x ~ 1 where 0 < x < 1:
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

    for j in ((n + 3)..i32::MAX).step_by(2) {
        p *= l2/(((j - 1)*j) as f64);
        let old_sum = sum;
        sum += zeta::zeta(n - j)*p;
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
    let mut sum = fac::fac(-n)*(-lnz).powi(n - 1);
    let (mut k, mut lnzk) = if is_even(n) {
        (1, lnz)
    } else {
        sum += zeta::zeta(n);
        (2, lnz2)
    };

    loop {
        let term = zeta::zeta(n - k)*fac::inv_fac(k)*lnzk;
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
fn li_series<T>(n: i32, z: T) -> T
    where T: AddAssign + Div<f64, Output = T> + Mul<Output = T> + MulAssign + PartialEq + Copy + IsInfinite<T>
{
    let mut sum = z;
    let mut zn = z*z;

    for k in 2..i32::MAX {
        let term = zn/(k as f64).powi(n);
        if term.is_inf() { break; }
        let old_sum = sum;
        sum += term;
        if sum == old_sum { break; }
        zn *= z;
    }

    sum
}

// Table[PolyLog[-2n+1,-1], {n,1,109}]
const LI_MINUS_1_COEFF_NEG: [f64; 109] = [
   -0.25, 0.125           , -0.25                  ,  1.0625                ,
   -7.75                  ,  86.375                , -1365.25               ,
    29049.03125           , -800572.75             ,  2.7741322625e7        ,
   -1.18052913025e9       ,  6.05239800516875e10   , -3.67941677853775e12   ,
    2.6170760990658388e014, -2.1531418140800295e016,  2.0288775575173016e018,
   -2.1708009902623771e020,  2.6173826968455815e022, -3.5324148876863878e024,
    5.3042033406864907e026, -8.8138218364311577e028,  1.6128065107490779e031,
   -3.2355470001722734e033,  7.0876727476537493e035, -1.6890450341293966e038,
    4.3639690731216831e040, -1.2185998827061261e043,  3.6670584803153006e045,
   -1.1859898526302099e048,  4.1120769493584015e050, -1.5249042436787620e053,
    6.0349693196941307e055, -2.5437161764210696e058,  1.1396923802632288e061,
   -5.4180861064753979e063,  2.7283654799994374e066, -1.4529750514918543e069,
    8.1705519371067450e071, -4.8445781606678368e074,  3.0246694206649519e077,
   -1.9858807961690493e080,  1.3694474620720087e083, -9.9070382984295808e085,
    7.5103780796592646e088, -5.9598418264260881e091,  4.9455988887500020e094,
   -4.2873596927020241e097,  3.8791952037716163e100, -3.6600317773156342e103,
    3.5978775704117284e106, -3.6818662617467813e109,  3.9192743066421374e112,
   -4.3363921885063858e115,  4.9833162711780838e118, -5.9438653020209606e121,
    7.3533439019770134e124, -9.4293465716973561e127,  1.2525196404154548e131,
   -1.7223787163994400e134,  2.4505178680729537e137, -3.6051616659014189e140,
    5.4813803836499771e143, -8.6083892012122616e146,  1.3957139354298160e150,
   -2.3350508860591630e153,  4.0291297374794860e156, -7.1669946227411534e159,
    1.3136385964069363e163, -2.4799083462304252e166,  4.8198083696385558e169,
   -9.6400031196958281e172,  1.9833611905147644e176, -4.1959717912682865e179,
    9.1243724595750010e182, -2.0386902382464212e186,  4.6786408066350383e189,
   -1.1024400389046488e193,  2.6662916424238258e196, -6.6165585014771755e199,
    1.6841726974970032e203, -4.3957474813006951e206,  1.1760766011899571e210,
   -3.2245094671360478e213,  9.0570855543185808e216, -2.6054618058433054e220,
    7.6741449421726560e223, -2.3136880427961752e227,  7.1382598572408242e230,
   -2.2530900128907084e234,  7.2736404696018159e237, -2.4010608416429639e241,
    8.1026279414941787e244, -2.7945745738098571e248,  9.8485095122481192e251,
   -3.5456055356238575e255,  1.3036999220919921e259, -4.8948166866453784e262,
    1.8761736309852136e266, -7.3399918877807488e269,  2.9303136033539038e273,
   -1.1935494277949469e277,  4.9589310621971370e280, -2.1012240064879845e284,
    9.0784179834777353e287, -3.9987113012775244e291,  1.7952380922182709e295,
   -8.2136799002055846e298,  3.8290431596908477e302, -1.8184610414701105e306
];

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
    if n < 0 {
        if is_even(n) {
            0.0
        } else if ((-(1 + n)/2) as usize) < LI_MINUS_1_COEFF_NEG.len() {
            LI_MINUS_1_COEFF_NEG[(-(1 + n)/2) as usize]
        } else if is_even((1 - n)/2) {
            std::f64::INFINITY
        } else {
            std::f64::NEG_INFINITY
        }
    } else if n == 0 {
        -0.5
    } else if n as usize <= LI_MINUS_1_COEFF.len() {
        LI_MINUS_1_COEFF[(n - 1) as usize]
    } else {
        -1.0
    }
}

#[test]
fn test_li_minus_1() {
    assert!(li_minus_1(-222) ==  0.0);
    assert!(li_minus_1(-221).is_infinite());
    assert!(li_minus_1(-220) ==  0.0);
    assert!(li_minus_1(-219).is_infinite());
    assert!(li_minus_1(-218) ==  0.0);
    assert!(li_minus_1(-217) == -1.8184610414701105e306);
    assert!(li_minus_1(  -4) ==  0.0);
    assert!(li_minus_1(  -3) ==  1.0/8.0);
    assert!(li_minus_1(  -2) ==  0.0);
    assert!(li_minus_1(  -1) == -0.25);
    assert!(li_minus_1(   1) == -0.69314718055994531);
    assert!(li_minus_1(   2) == -0.82246703342411322);
    assert!(li_minus_1(  52) == -0.9999999999999998);
    assert!(li_minus_1(  53) == -0.9999999999999999);
    assert!(li_minus_1(  54) == -0.9999999999999999);
    assert!(li_minus_1(  55) == -1.0);
    assert!(li_minus_1(  56) == -1.0);
}

/// provides the is_inf() function to test whether T is non-finite
/// (nan or infinite)
trait IsInfinite<T> {
    fn is_inf(&self) -> bool;
}

impl IsInfinite<f64> for f64 {
    fn is_inf(&self) -> bool {
        !self.is_finite()
    }
}

impl IsInfinite<Complex<f64>> for Complex<f64> {
    fn is_inf(&self) -> bool {
        !self.is_finite()
    }
}
