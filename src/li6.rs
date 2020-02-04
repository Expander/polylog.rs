use std;
use num::complex::Complex;

/// Provides the sixths order polylogarithm function `li6()` of a
/// number of type `T`.
pub trait Li6<T> {
    fn li6(&self) -> T;
}

impl Li6<Complex<f64>> for Complex<f64> {
    /// Returns the sixths order polylogarithm of a complex number of type
    /// `Complex<f64>`.
    ///
    /// # Example:
    /// ```
    /// extern crate num;
    /// extern crate polylog;
    /// use num::complex::Complex;
    /// use polylog::Li6;
    ///
    /// fn main() {
    ///     let z = Complex::new(1.0, 1.0);
    ///     println!("Li6({}) = {}", z, z.li6());
    /// }
    /// ```
    fn li6(&self) -> Complex<f64> {
        let pi  = 3.141592653589793;
        let pi2 = pi*pi;
        let eps = std::f64::EPSILON;
        let z6  = 1.017343061984449; // zeta(6)
        let bf  = [
            1.                   , -31./64.              ,
            1.524134087791495e-01, -3.436555587705761e-02,
            5.717479723936899e-03, -6.818045374657064e-04,
            4.996036194873449e-05, -4.916605119603904e-07,
           -3.063297516130216e-07,  1.441459927084909e-08,
            3.727243823092410e-09, -3.730086734548760e-10,
           -5.124652681608583e-11,  9.054193095663668e-12,
            6.738188261551251e-13, -2.121583115030313e-13,
           -6.840881171901169e-15,  4.869117846200558e-15
        ];

        if is_close(self, 0.0, eps) {
            return Complex::new(0.0, 0.0);
        }
        if is_close(self, 1.0, eps) {
            return Complex::new(z6, 0.0);
        }
        if is_close(self, -1.0, eps) {
            return Complex::new(-31./32.*z6, 0.0);
        }

        let az  = self.norm();
        let pz  = self.arg();
        let lnz = az.ln();

        if lnz*lnz + pz*pz < 1.0 { // |log(z)| < 1
            let u  = Complex::new(lnz, pz);
            let u2 = u*u;
            let c1 = 1.036927755143370; // zeta(5)
            let c2 = 0.5411616168555691;
            let c3 = 0.2003428171932657;
            let c4 = 0.06853891945200943;
            let c5 = (137./60. - (-u).cln())/120.;
            let c6 = -1./1440.;

            let cs = [
                -1.653439153439153e-05, 2.296443268665491e-08,
                -9.941312851365761e-11, 6.691268265342339e-13,
                -5.793305857439255e-15
            ];

            return z6 + u * c1 +
                u2 * (c2 + u * c3 +
                u2 * (c4 + u * c5 +
                u2 * (c6 +
                u * (cs[0] +
                u2 * (cs[1] +
                u2 * (cs[2] +
                u2 * (cs[3] +
                u2 * (cs[4]))))))));
        }

        let (u, rest, sgn) = if az <= 1. {
            (-(1.0 - self).cln(), Complex::new(0.0, 0.0), 1.)
        } else { // az > 1.
            let pi4 = pi2*pi2;
            let pi6 = pi2*pi4;
            let arg = {
                let res = pi + pz;
                if res > pi { res - 2.*pi } else { res }
            };
            let lmz = Complex::new(lnz, arg); // (-self).cln()
            let lmz2 = pow2(lmz);
            (-(1. - 1./self).cln(), -31.*pi6/15120. + lmz2*(-7./720.*pi4 + lmz2*(-1./144.*pi2 - 1./720.*lmz2)), -1.)
        };

        rest + sgn * (
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
            u * (bf[17])))))))))))))))))))
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
