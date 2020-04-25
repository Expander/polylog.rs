use std;
use num::complex::Complex;

/// Provides the trilogarithm function `li3()` of a number of type
/// `T`.
pub trait Li3<T> {
    fn li3(&self) -> T;
}

impl Li3<Complex<f64>> for Complex<f64> {
    /// Returns the trilogarithm of a complex number of type `Complex<f64>`.
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
    fn li3(&self) -> Complex<f64> {
        let pi  = 3.141592653589793;
        let pi2 = pi*pi;
        let eps = std::f64::EPSILON;
        let z2  = 1.644934066848226;
        let z3  = 1.202056903159594;
        let bf  = [
            1., -3./8., 17./216., -5./576.,
            1.296296296296296e-04,  8.101851851851851e-05,
           -3.419357160853759e-06, -1.328656462585034e-06,
            8.660871756109851e-08,  2.526087595532039e-08,
           -2.144694468364064e-09, -5.140110622012978e-10,
            5.249582114600829e-11,  1.088775440663631e-11,
           -1.277939609449369e-12, -2.369824177308745e-13,
            3.104357887965462e-14,  5.261758629912506e-15,
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

        let nz  = self.norm_sqr();
        let pz  = self.arg();
        let lnz = 0.5*nz.ln();

        if lnz*lnz + pz*pz < 1. { // |log(z)| < 1
            let u  = Complex::new(lnz, pz);
            let u2 = u*u;
            let c0 = z3 + u*(z2 - u2/12.);
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

        let (u, rest) = if nz <= 1.0 {
            (-(1. - self).cln(), Complex::new(0.,0.))
        } else { // nz > 1.0
            let arg = if pz > 0.0 { pz - pi } else { pz + pi };
            let lmz = Complex::new(lnz, arg); // (-self).cln()
            (-(1. - 1./self).cln(), -lmz*(sqr(lmz)/6. + pi2/6.))
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

fn sqr(z : Complex<f64>) -> Complex<f64> {
    z * z
}

trait CLn<T> {
    fn cln(&self) -> T;
}

impl CLn<Complex<f64>> for Complex<f64> {
    fn cln(&self) -> Complex<f64> {
        let z = Complex::new(
            if self.re == 0. { 0. } else { self.re },
            if self.im == 0. { 0. } else { self.im },
        );
        Complex::new(0.5*z.norm_sqr().ln(), z.arg())
    }
}
