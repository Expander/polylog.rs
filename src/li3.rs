use num::complex::Complex;
use crate::cln::CLn;

/// Provides the 3rd order polylogarithm (trilogarithm) function
/// `li3()` of a number of type `T`.
pub trait Li3<T> {
    fn li3(&self) -> T;
}

impl Li3<f64> for f64 {
    /// Returns the real trilogarithm of a real number of type `f64`.
    ///
    /// Implemented as rational function approximations with a maximum
    /// error of less than 2.050e-17 [[arXiv:2308.11619]].
    ///
    /// [arXiv:2308.11619]: https://arxiv.org/abs/2308.11619
    ///
    /// # Example:
    /// ```
    /// use polylog::Li3;
    ///
    /// assert!((1.0_f64.li3() - 1.2020569031595943_f64).abs() < std::f64::EPSILON);
    /// ```
    fn li3(&self) -> f64 {
        let z2 = 1.6449340668482264;
        let z3 = 1.2020569031595943;
        let x = *self;

        // transformation to [-1,0] and [0,1/2]
        if x < -1.0 {
            let l = (-x).ln();
            li3_neg(1.0/x) - l*(z2 + 1.0/6.0*l*l)
        } else if x == -1.0 {
            -0.75*z3
        } else if x < 0.0 {
            li3_neg(x)
        } else if x == 0.0 {
            x
        } else if x < 0.5 {
            li3_pos(x)
        } else if x == 0.5 {
            0.53721319360804020
        } else if x < 1.0 {
            let l = x.ln();
            -li3_neg(1.0 - x.recip()) - li3_pos(1.0 - x) + z3 + l*(z2 + l*(-0.5*(-x).ln_1p() + 1.0/6.0*l))
        } else if x == 1.0 {
            z3
        } else if x < 2.0 {
            let l = x.ln();
            -li3_neg(1.0 - x) - li3_pos(1.0 - x.recip()) + z3 + l*(z2 + l*(-0.5*(x - 1.0).ln() + 1.0/6.0*l))
        } else { // x >= 2.0
            let l = x.ln();
            li3_pos(x.recip()) + l*(2.0*z2 - 1.0/6.0*l*l)
        }
    }
}

// Li_3(x) for x in [-1,0]
fn li3_neg(x: f64) -> f64 {
    let cp = [
        0.9999999999999999795e+0, -2.0281801754117129576e+0,
        1.4364029887561718540e+0, -4.2240680435713030268e-1,
        4.7296746450884096877e-2, -1.3453536579918419568e-3
    ];
    let cq = [
        1.0000000000000000000e+0, -2.1531801754117049035e+0,
        1.6685134736461140517e+0, -5.6684857464584544310e-1,
        8.1999463370623961084e-2, -4.0756048502924149389e-3,
        3.4316398489103212699e-5
    ];

    let x2 = x*x;
    let x4 = x2*x2;
    let p = cp[0] + x * cp[1] + x2 * (cp[2] + x * cp[3]) +
            x4 * (cp[4] + x * cp[5]);
    let q = cq[0] + x * cq[1] + x2 * (cq[2] + x * cq[3]) +
            x4 * (cq[4] + x * cq[5] + x2 * cq[6]);

    x*p/q
}

// Li_3(x) for x in [0,1/2]
fn li3_pos(x: f64) -> f64 {
    let cp = [
        0.9999999999999999893e+0, -2.5224717303769789628e+0,
        2.3204919140887894133e+0, -9.3980973288965037869e-1,
        1.5728950200990509052e-1, -7.5485193983677071129e-3
    ];
    let cq = [
        1.0000000000000000000e+0, -2.6474717303769836244e+0,
        2.6143888433492184741e+0, -1.1841788297857667038e+0,
        2.4184938524793651120e-1, -1.8220900115898156346e-2,
        2.4927971540017376759e-4
    ];

    let x2 = x*x;
    let x4 = x2*x2;
    let p = cp[0] + x * cp[1] + x2 * (cp[2] + x * cp[3]) +
            x4 * (cp[4] + x * cp[5]);
    let q = cq[0] + x * cq[1] + x2 * (cq[2] + x * cq[3]) +
            x4 * (cq[4] + x * cq[5] + x2 * cq[6]);

    x*p/q
}

impl Li3<Complex<f64>> for Complex<f64> {
    /// Returns the trilogarithm of a complex number of type `Complex<f64>`.
    ///
    /// # Example:
    /// ```
    /// use num::complex::Complex;
    /// use polylog::Li3;
    ///
    /// assert!((Complex::new(1.0_f64, 1.0_f64).li3() - Complex::new(0.8711588834109380_f64, 1.2670834418889240_f64)).norm() < 2.0_f64*std::f64::EPSILON);
    /// ```
    fn li3(&self) -> Complex<f64> {
        let pi  = std::f64::consts::PI;
        let z2  = 1.6449340668482264;
        let z3  = 1.2020569031595943;

        if self.im == 0.0 {
            if self.re <= 1.0 {
                Complex::new(self.re.li3(), self.im)
            } else { // rz > 1.0
                let l = self.re.ln();
                Complex::new(self.re.li3(), -0.5*pi*l*l)
            }
        } else {
            let nz  = self.norm();
            let pz  = self.arg();
            let lnz = nz.ln();

            if lnz*lnz + pz*pz < 1.0 { // |log(z)| < 1
                let u  = Complex::new(lnz, pz);
                let u2 = u*u;
                let u4 = u2*u2;
                let u8 = u4*u4;
                let c0 = z3 + u*(z2 - u2/12.0);
                let c1 = 0.25 * (3.0 - 2.0*(-u).cln());

                let cs = [
                    -3.4722222222222222e-03, 1.1574074074074074e-05,
                    -9.8418997228521038e-08, 1.1482216343327454e-09,
                    -1.5815724990809166e-11, 2.4195009792525152e-13,
                    -3.9828977769894877e-15
                ];

                c0 +
                c1*u2 +
                u4*(cs[0] + u2*cs[1]) +
                u8*(cs[2] + u2*cs[3] + u4*(cs[4] + u2*cs[5])) +
                u8*u8*cs[6]
            } else if nz <= 1.0 {
                cli3_unit_circle(-(1.0 - self).cln())
            } else { // nz > 1
                let arg = if pz > 0.0 { pz - pi } else { pz + pi };
                let lmz = Complex::new(lnz, arg); // (-self).cln()
                cli3_unit_circle(-(1.0 - 1.0/self).cln()) - lmz*(lmz*lmz/6.0 + z2)
            }
        }
    }
}

/// series approximation of Li3(z) for |z| <= 1
/// in terms of x = -ln(1 - z)
fn cli3_unit_circle(x: Complex<f64>) -> Complex<f64> {
    let bf  = [
        1.0, -3.0/8.0, 17.0/216.0, -5.0/576.0,
        1.2962962962962963e-04,  8.1018518518518519e-05,
       -3.4193571608537595e-06, -1.3286564625850340e-06,
        8.6608717561098513e-08,  2.5260875955320400e-08,
       -2.1446944683640648e-09, -5.1401106220129789e-10,
        5.2495821146008294e-11,  1.0887754406636318e-11,
       -1.2779396094493695e-12, -2.3698241773087452e-13,
        3.1043578879654623e-14,  5.2617586299125061e-15,
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
