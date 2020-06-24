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
        let pi  = std::f64::consts::PI;
        let pi2 = pi*pi;
        let z5  = 1.0369277551433699; // zeta(5)
        let bf  = [
            1.                    , -15./32.               ,
            1.3953189300411523e-01, -2.8633777006172840e-02,
            4.0317412551440329e-03, -3.3985018004115226e-04,
            4.5445184621617666e-06,  2.3916808048569012e-06,
           -1.2762692600122747e-07, -3.1628984306505932e-08,
            3.2848118445335192e-09,  4.7613713995660579e-10,
           -8.0846898171909830e-11, -7.2387648587737207e-12,
            1.9439760115173968e-12,  1.0256978405977236e-13,
           -4.6180551009884830e-14, -1.1535857196470580e-15,
            1.0903545401333394e-15
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
            let c1 = 1.0823232337111382; // zeta(4)
            let c2 = 0.60102845157979714; // zeta(3)/2
            let c3 = 0.27415567780803774;
            let c4 = (25.0/12.0 - (-u).cln())/24.0;
            let c5 = -1.0/240.0;

            let cs = [
                -1.1574074074074074e-04, 2.0667989417989418e-07,
                -1.0935444136502338e-09, 8.6986487449450412e-12,
                -8.6899587861588824e-14, 1.0081254080218813e-15
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

        let u2 = u*u;
        let u4 = u2*u2;
        let u8 = u4*u4;

        rest +
        u*bf[0] +
        u2*(bf[1] + u*bf[2]) +
        u4*(bf[3] + u*bf[4] + u2*(bf[5] + u*bf[6])) +
        u8*(bf[7] + u*bf[8] + u2*(bf[9] + u*bf[10]) +
            u4*(bf[11] + u*bf[12] + u2*(bf[13] + u*bf[14]))) +
        u8*u8*(bf[15] + u*bf[16] + u2*(bf[17] + u*bf[18]))
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
