use num::complex::Complex;

/// Provides the fifth order polylogarithm function `li5()` of a
/// number of type `T`.
pub trait Li5<T> {
    fn li5(&self) -> T;
}

impl Li5<Complex<f64>> for Complex<f64> {
    /// Returns the fifth order polylogarithm of a complex number of type
    /// `Complex<f64>`.
    ///
    /// # Example:
    /// ```
    /// extern crate num;
    /// extern crate polylog;
    /// use num::complex::Complex;
    /// use polylog::Li5;
    ///
    /// fn main() {
    ///     let z = Complex::new(1.0, 1.0);
    ///     println!("Li5({}) = {}", z, z.li5());
    /// }
    /// ```
    fn li5(&self) -> Complex<f64> {
        let pi  = 3.141592653589793;
        let pi2 = pi*pi;
        let z5  = 1.036927755143370; // zeta(5)
        let bf  = [
            1.                   , -15./32.              ,
            1.395318930041152e-01, -2.863377700617283e-02,
            4.031741255144032e-03, -3.398501800411522e-04,
            4.544518462161766e-06,  2.391680804856901e-06,
           -1.276269260012274e-07, -3.162898430650593e-08,
            3.284811844533519e-09,  4.761371399566057e-10,
           -8.084689817190983e-11, -7.238764858773720e-12,
            1.943976011517396e-12,  1.025697840597723e-13,
           -4.618055100988483e-14, -1.153585719647058e-15,
            1.090354540133339e-15
        ];

        if self.im == 0.0 {
            if self.re == 0.0 {
                return Complex::new(0., 0.);
            }
            if self.re == 1.0 {
                return Complex::new(z5, 0.);
            }
            if self.re == -1.0 {
                return Complex::new(-15./16.*z5, 0.0);
            }
        }

        let nz  = self.norm_sqr();
        let pz  = self.arg();
        let lnz = 0.5*nz.ln();

        if lnz*lnz + pz*pz < 1.0 { // |log(z)| < 1
            let u  = Complex::new(lnz, pz);
            let u2 = u*u;
            let c1 = 1.082323233711138; // zeta(4)
            let c2 = 0.6010284515797971; // zeta(3)/2
            let c3 = 0.2741556778080377;
            let c4 = (25.0/12.0 - (-u).cln())/24.0;
            let c5 = -1.0/240.0;

            let cs = [
                -1.157407407407407e-04, 2.066798941798942e-07,
                -1.093544413650234e-09, 8.698648744945041e-12,
                -8.689958786158882e-14, 1.008125408021881e-15
            ];

            return z5 + u * c1 +
                u2 * (c2 + u * c3 +
                u2 * (c4 + u * c5 +
                u2 * (cs[0] +
                u2 * (cs[1] +
                u2 * (cs[2] +
                u2 * (cs[3] +
                u2 * (cs[4] +
                u2 * (cs[5]))))))));
        }

        let (u, rest) = if nz <= 1.0 {
            (-(1.0 - self).cln(), Complex::new(0.0, 0.0))
        } else { // nz > 1.0
            let pi4  = pi2*pi2;
            let arg = if pz > 0.0 { pz - pi } else { pz + pi };
            let lmz = Complex::new(lnz, arg); // (-self).cln()
            let lmz2 = lmz*lmz;
            (-(1. - 1./self).cln(), -1./360.*lmz*(7.*pi4 + lmz2*(10.*pi2 + 3.*lmz2)))
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
            u * (bf[17] +
            u * (bf[18])))))))))))))))))))
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
