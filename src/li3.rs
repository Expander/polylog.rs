use std;
use num::complex::Complex;

/// Provides the trilogarithm function `li3()` of a number of type
/// `T`.
pub trait Li3<T> {
    fn li3(&self) -> T;
}

/// Returns the trilogarithm of a real number of type `f64`.
///
/// # Example:
/// ```
/// extern crate num;
/// extern crate polylog;
/// use num::complex::Complex;
/// use polylog::Li3;
///
/// fn main() {
///     let z = Complex::new(1.0, 1.0);
///     println!("Li3({}) = {}", z, z.li3());
/// }
/// ```
impl Li3<Complex<f64>> for Complex<f64> {
    fn li3(&self) -> Complex<f64> {
        let pi  = 3.141592653589793;
        let pi2 = pi*pi;
        let eps = std::f64::EPSILON;
        let z2  = 1.644934066848226;
        let z3  = 1.202056903159594;
        let bf  = [
            1., -3./8., 17./216., -5./576.,
             0.00012962962962962962962962962962963,  0.000081018518518518518518518518518519,
            -3.4193571608537594932152755282007e-06, -1.3286564625850340136054421768707e-06 ,
             8.6608717561098513479465860418241e-08,  2.5260875955320399764844209288654e-08 ,
            -2.1446944683640647609338850757365e-09, -5.1401106220129789153358176927200e-10 ,
             5.2495821146008294363940888085581e-11,  1.0887754406636318375372971570425e-11 ,
            -1.2779396094493695305581831754072e-12, -2.3698241773087452099797778810124e-13 ,
             3.1043578879654622942847532704656e-14,  5.2617586299125060841318392511225e-15 ,
        ];

        if is_close(self, 0., eps) {
            return Complex::new(0., 0.);
        }
        if is_close(self, 1., eps) {
            return Complex::new(z3, 0.);
        }
        if is_close(self, -1., eps) {
            return Complex::new(-0.75*z3, 0.);
        }
        if is_close(self, 0.5, eps) {
            let ln2  = 0.6931471805599453; // ln(2)
            let ln23 = 0.3330246519889295; // ln(2)^3
            return Complex::new((-2.*pi2*ln2 + 4.*ln23 + 21.*z3)/24., 0.);
        }

        let az  = self.norm();
        let pz  = self.arg();
        let lnz = az.ln();

        if lnz*lnz + pz*pz < 1. { // |log(z)| < 1
           let u  = self.cln();
           let u2 = u*u;
           let u3 = u*u2;
           let c0 = z3 + z2*u - u3/12.;
           let c1 = 0.25 * (3.0 - 2.0*(-u).cln());

           let cs = [
              -3.472222222222222e-03, 1.157407407407407e-05,
              -9.841899722852104e-08, 1.148221634332745e-09,
              -1.581572499080917e-11, 2.419500979252515e-13,
              -3.982897776989488e-15
           ];

           return c0 +
              u2 * (c1 +
              u2 * (cs[0] +
              u2 * (cs[1] +
              u2 * (cs[2] +
              u2 * (cs[3] +
              u2 * (cs[4] +
              u2 * (cs[5] +
              u2 * (cs[6]))))))));
        }

        let (u, rest) = if az <= 1. {
            (-(1. - self).cln(), Complex::new(0.,0.))
        } else { // az > 1.
            let lmz = (-self).cln();
            (-(1. - 1./self).cln(), -lmz*(pow2(lmz)/6. + pi2/6.))
        };

        rest +
        u * (bf[0] +
        u * (bf[1] +
        u * (bf[2] +
        u * (bf[3] +
        u * (bf[4] +
        u * (bf[5] +
        u * (bf[6] +
        u * (bf[7] +
        u * (bf[8] +
        u * (bf[9] +
        u * (bf[10] +
        u * (bf[11] +
        u * (bf[12] +
        u * (bf[13] +
        u * (bf[14] +
        u * (bf[15] +
        u * (bf[16] +
        u * (bf[17]))))))))))))))))))
    }
}

fn is_close(a : &Complex<f64>, b : f64, eps : f64) -> bool {
    (a.re - b).abs() < eps && (a.im).abs() < eps
}

fn pow2(z : Complex<f64>) -> Complex<f64> {
    z * z
}

trait CLn<T> {
    fn cln(&self) -> T;
}

impl CLn<Complex<f64>> for Complex<f64> {
    fn cln(&self) -> Complex<f64> {
        Complex::new(
            if self.re == 0. { 0. } else { self.re },
            if self.im == 0. { 0. } else { self.im },
        ).ln()
    }
}
