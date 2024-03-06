use num::complex::Complex;
use crate::cln::CLn;

/// Provides the 5-th order polylogarithm function `li5()` of a
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
    /// use num::complex::Complex;
    /// use polylog::Li5;
    ///
    /// assert!((Complex::new(1.0_f64, 1.0_f64).li5() - Complex::new(0.9874666591701124_f64, 1.0684416071074221_f64)).norm() < 2.0_f64*std::f64::EPSILON);
    /// ```
    fn li5(&self) -> Complex<f64> {
        let pi  = std::f64::consts::PI;
        let pi2 = pi*pi;
        let z5  = 1.0369277551433699; // zeta(5)

        if self.im == 0.0 && self.re == 0.0 {
            *self
        } else if self.im == 0.0 && self.re == 1.0 {
            Complex::new(z5, self.im)
        } else if self.im == 0.0 && self.re == -1.0 {
            Complex::new(-15.0/16.0*z5, self.im)
        } else {
            let nz  = self.norm();
            let pz  = self.arg();
            let lnz = nz.ln();

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

                z5 + u * c1 +
                u2 * (c2 + u * c3 +
                u2 * (c4 + u * c5 +
                u2 * (cs[0] +
                u2 * (cs[1] +
                u2 * (cs[2] +
                u2 * (cs[3] +
                u2 * (cs[4] +
                u2 * (cs[5]))))))))
            } else if nz <= 1.0 {
                cli5_unit_circle(-(1.0 - self).cln())
            } else { // nz > 1.0
                let pi4  = pi2*pi2;
                let arg = if pz > 0.0 { pz - pi } else { pz + pi };
                let lmz = Complex::new(lnz, arg); // (-self).cln()
                let lmz2 = lmz*lmz;
                cli5_unit_circle(-(1.0 - 1.0/self).cln()) - 1.0/360.0*lmz*(7.0*pi4 + lmz2*(10.0*pi2 + 3.0*lmz2))
            }
        }
    }
}

/// series approximation of Li5(z) for |z| <= 1
/// in terms of x = -ln(1 - z)
fn cli5_unit_circle(x: Complex<f64>) -> Complex<f64> {
    let bf  = [
        1.0                   , -15.0/32.0             ,
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

    let x2 = x*x;
    let x4 = x2*x2;
    let x8 = x4*x4;

    x*bf[0] +
    x2*(bf[1] + x*bf[2]) +
    x4*(bf[3] + x*bf[4] + x2*(bf[5] + x*bf[6])) +
    x8*(bf[7] + x*bf[8] + x2*(bf[9] + x*bf[10]) +
        x4*(bf[11] + x*bf[12] + x2*(bf[13] + x*bf[14]))) +
    x8*x8*(bf[15] + x*bf[16] + x2*(bf[17] + x*bf[18]))
}
