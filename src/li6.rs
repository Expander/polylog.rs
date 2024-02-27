use num::complex::Complex;
use crate::cln::CLn;

/// Provides the 6-th order polylogarithm function `li6()` of a
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
    /// use num::complex::Complex;
    /// use polylog::Li6;
    ///
    /// assert!((Complex::new(1.0_f64, 1.0_f64).li6() - Complex::new(0.9961497968353170_f64, 1.0335544477237482_f64)).norm() < 2.0_f64*std::f64::EPSILON);
    /// ```
    fn li6(&self) -> Complex<f64> {
        let pi  = std::f64::consts::PI;
        let pi2 = pi*pi;
        let z6  = 1.0173430619844491; // zeta(6)

        if self.im == 0.0 && self.re == 0.0 {
            Complex::new(0.0, 0.0)
        } else if self.im == 0.0 && self.re == 1.0 {
            Complex::new(z6, 0.0)
        } else if self.im == 0.0 && self.re == -1.0 {
            Complex::new(-31.0/32.0*z6, 0.0)
        } else {
            let nz  = self.norm();
            let pz  = self.arg();
            let lnz = nz.ln();

            if lnz*lnz + pz*pz < 1.0 { // |log(z)| < 1
                let u  = Complex::new(lnz, pz);
                let u2 = u*u;
                let c1 = 1.0369277551433699; // zeta(5)
                let c2 = 0.54116161685556910;
                let c3 = 0.20034281719326571;
                let c4 = 0.068538919452009435;
                let c5 = (137.0/60.0 - (-u).cln())/120.0;
                let c6 = -1.0/1440.0;

                let cs = [
                    -1.6534391534391534e-05, 2.2964432686654909e-08,
                    -9.9413128513657614e-11, 6.6912682653423394e-13,
                    -5.7933058574392549e-15
                ];

                z6 + u * c1 +
                u2 * (c2 + u * c3 +
                u2 * (c4 + u * c5 +
                u2 * (c6 +
                u * (cs[0] +
                u2 * (cs[1] +
                u2 * (cs[2] +
                u2 * (cs[3] +
                u2 * (cs[4]))))))))
            } else if nz <= 1.0 {
                cli6_unit_circle(-(-self).cln_1p())
            } else { // nz > 1.0
                let pi4 = pi2*pi2;
                let pi6 = pi2*pi4;
                let arg = if pz > 0.0 { pz - pi } else { pz + pi };
                let lmz = Complex::new(lnz, arg); // (-self).cln()
                let lmz2 = lmz*lmz;
                -cli6_unit_circle(-(-1.0/self).cln_1p()) - 31.0*pi6/15120.0 + lmz2*(-7.0/720.0*pi4 + lmz2*(-1.0/144.0*pi2 - 1.0/720.0*lmz2))
            }
        }
    }
}

/// series approximation of Li6(z) for |z| <= 1
/// in terms of x = -ln(1 - z)
fn cli6_unit_circle(x: Complex<f64>) -> Complex<f64> {
    let bf  = [
        1.0                   , -31.0/64.0             ,
        1.5241340877914952e-01, -3.4365555877057613e-02,
        5.7174797239368999e-03, -6.8180453746570645e-04,
        4.9960361948734493e-05, -4.9166051196039048e-07,
       -3.0632975161302164e-07,  1.4414599270849095e-08,
        3.7272438230924107e-09, -3.7300867345487607e-10,
       -5.1246526816085832e-11,  9.0541930956636683e-12,
        6.7381882615512517e-13, -2.1215831150303135e-13,
       -6.8408811719011698e-15,  4.8691178462005581e-15
    ];

    let x2 = x*x;
    let x4 = x2*x2;
    let x8 = x4*x4;

    x*bf[0] +
    x2*(bf[1] + x*bf[2]) +
    x4*(bf[3] + x*bf[4] + x2*(bf[5] + x*bf[6])) +
    x8*(bf[7] + x*bf[8] + x2*(bf[9] + x*bf[10]) +
        x4*(bf[11] + x*bf[12] + x2*(bf[13] + x*bf[14]))) +
    x8*x8*(bf[15] + x*bf[16] + x2*bf[17])
}
