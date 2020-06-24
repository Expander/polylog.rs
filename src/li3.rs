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
        let pi  = std::f64::consts::PI;
        let z2  = 1.6449340668482264;
        let z3  = 1.2020569031595943;
        let bf  = [
            1., -3./8., 17./216., -5./576.,
            1.2962962962962963e-04,  8.1018518518518519e-05,
           -3.4193571608537595e-06, -1.3286564625850340e-06,
            8.6608717561098513e-08,  2.5260875955320400e-08,
           -2.1446944683640648e-09, -5.1401106220129789e-10,
            5.2495821146008294e-11,  1.0887754406636318e-11,
           -1.2779396094493695e-12, -2.3698241773087452e-13,
            3.1043578879654623e-14,  5.2617586299125061e-15,
        ];

        if self.im == 0.0 {
            if self.re == 0.0 {
                return Complex::new(0., 0.);
            }
            if self.re == 1.0 {
                return Complex::new(z3, 0.);
            }
            if self.re == -1.0 {
                return Complex::new(-0.75*z3, 0.);
            }
            if self.re == 0.5 {
                return Complex::new(0.53721319360804020, 0.);
            }
        }

        let nz  = self.norm_sqr();
        let pz  = self.arg();
        let lnz = 0.5*nz.ln();

        if lnz*lnz + pz*pz < 1. { // |log(z)| < 1
            let u  = Complex::new(lnz, pz);
            let u2 = u*u;
            let u4 = u2*u2;
            let u8 = u4*u4;
            let c0 = z3 + u*(z2 - u2/12.);
            let c1 = 0.25 * (3.0 - 2.0*(-u).cln());

            let cs = [
                -3.4722222222222222e-03, 1.1574074074074074e-05,
                -9.8418997228521038e-08, 1.1482216343327454e-09,
                -1.5815724990809166e-11, 2.4195009792525152e-13,
                -3.9828977769894877e-15
            ];

            return c0 +
               c1*u2 +
               u4*(cs[0] + u2*cs[1]) +
               u8*(cs[2] + u2*cs[3] + u4*(cs[4] + u2*cs[5])) +
               u8*u8*cs[6];
        }

        let (u, rest) = if nz <= 1.0 {
            (-(1. - self).cln(), Complex::new(0.,0.))
        } else { // nz > 1.0
            let arg = if pz > 0.0 { pz - pi } else { pz + pi };
            let lmz = Complex::new(lnz, arg); // (-self).cln()
            (-(1. - 1./self).cln(), -lmz*(lmz*lmz/6. + z2))
        };

        let u2 = u*u;
        let u4 = u2*u2;
        let u8 = u4*u4;

        rest +
        u*bf[0] +
        u2*(bf[1] + u*bf[2]) +
        u4*(bf[3] + u*bf[4] + u2*(bf[5] + u*bf[6])) +
        u8*(bf[7] + u*bf[8] + u2*(bf[9] + u*bf[10]) +
            u4*(bf[11] + u*bf[12] + u2*(bf[13] + u*bf[14]))) +
        u8*u8*(bf[15] + u*bf[16] + u2*bf[17])
    }
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
