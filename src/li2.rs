use num::complex::Complex;
use crate::cln::CLn;

/// Provides the 2nd order polylogarithm (dilogarithm) function
/// `li2()` of a number of type `T`.
pub trait Li2<T> {
    fn li2(&self) -> T;
}

/// rational function approximation of Re[Li2(x)] for x in [0, 1/2]
fn li2_approx_f32(x: f32) -> f32 {
    let cp = [ 1.00000020_f32, -0.780790946_f32, 0.0648256871_f32 ];
    let cq = [ 1.00000000_f32, -1.03077545_f32, 0.211216710_f32 ];

    let p = cp[0] + x*(cp[1] + x*cp[2]);
    let q = cq[0] + x*(cq[1] + x*cq[2]);

    x*p/q
}

/// rational function approximation of Re[Li2(x)] for x in [0, 1/2]
fn li2_approx_f64(x: f64) -> f64 {
    let cp = [
         0.9999999999999999502e+0,
        -2.6883926818565423430e+0,
         2.6477222699473109692e+0,
        -1.1538559607887416355e+0,
         2.0886077795020607837e-1,
        -1.0859777134152463084e-2
    ];
    let cq = [
         1.0000000000000000000e+0,
        -2.9383926818565635485e+0,
         3.2712093293018635389e+0,
        -1.7076702173954289421e+0,
         4.1596017228400603836e-1,
        -3.9801343754084482956e-2,
         8.2743668974466659035e-4
    ];

    let x2 = x*x;
    let x4 = x2*x2;
    let p = cp[0] + x*cp[1] + x2*(cp[2] + x*cp[3]) +
        x4*(cp[4] + x*cp[5]);
    let q = cq[0] + x*cq[1] + x2*(cq[2] + x*cq[3]) +
        x4*(cq[4] + x*cq[5] + x2*cq[6]);

    x*p/q
}

impl Li2<f32> for f32 {
    /// Returns the real dilogarithm of a real number of type `f32`.
    ///
    /// Implemented as rational function approximation with a maximum
    /// error of 5e-17 [[arXiv:2201.01678]].
    ///
    /// [arXiv:2201.01678]: https://arxiv.org/abs/2201.01678
    ///
    /// # Example:
    /// ```
    /// use polylog::Li2;
    ///
    /// assert!((1.0_f32.li2() - 1.64493407_f32).abs() < 2.0_f32*std::f32::EPSILON);
    /// ```
    fn li2(&self) -> f32 {
        let pi = std::f32::consts::PI;
        let x = *self;

        // transform to [0, 1/2]
        if x < -1.0_f32 {
            let l = (1.0_f32 - x).ln();
            li2_approx_f32(1.0_f32/(1.0_f32 - x)) - pi*pi/6.0_f32 + l*(0.5_f32*l - (-x).ln())
        } else if x == -1.0_f32 {
            -pi*pi/12.0_f32
        } else if x < 0.0_f32 {
            let l = (-x).ln_1p();
            -li2_approx_f32(x/(x - 1.0_f32)) - 0.5_f32*l*l
        } else if x == 0.0_f32 {
            0.0_f32
        } else if x < 0.5_f32 {
            li2_approx_f32(x)
        } else if x < 1.0_f32 {
            -li2_approx_f32(1.0_f32 - x) + pi*pi/6.0_f32 - x.ln()*(-x).ln_1p()
        } else if x == 1.0_f32 {
            pi*pi/6.0_f32
        } else if x < 2.0_f32 {
            let l = x.ln();
            li2_approx_f32(1.0_f32 - 1.0_f32/x) + pi*pi/6.0_f32 - l*((1.0_f32 - 1.0_f32/x).ln() + 0.5_f32*l)
        } else {
            let l = x.ln();
            -li2_approx_f32(1.0_f32/x) + pi*pi/3.0_f32 - 0.5_f32*l*l
        }
    }
}

impl Li2<f64> for f64 {
    /// Returns the real dilogarithm of a real number of type `f64`.
    ///
    /// Implemented as rational function approximation with a maximum
    /// error of 5e-17 [[arXiv:2201.01678]].
    ///
    /// [arXiv:2201.01678]: https://arxiv.org/abs/2201.01678
    ///
    /// # Example:
    /// ```
    /// use polylog::Li2;
    ///
    /// assert!((1.0_f64.li2() - 1.6449340668482264_f64).abs() < 2.0_f64*std::f64::EPSILON);
    /// ```
    fn li2(&self) -> f64 {
        let pi = std::f64::consts::PI;
        let x = *self;

        // transform to [0, 1/2]
        if x < -1.0_f64 {
            let l = (1.0_f64 - x).ln();
            li2_approx_f64(1.0_f64/(1.0_f64 - x)) - pi*pi/6.0_f64 + l*(0.5_f64*l - (-x).ln())
        } else if x == -1.0_f64 {
            -pi*pi/12.0_f64
        } else if x < 0.0_f64 {
            let l = (-x).ln_1p();
            -li2_approx_f64(x/(x - 1.0_f64)) - 0.5_f64*l*l
        } else if x == 0.0_f64 {
            0.0_f64
        } else if x < 0.5_f64 {
            li2_approx_f64(x)
        } else if x < 1.0_f64 {
            -li2_approx_f64(1.0_f64 - x) + pi*pi/6.0_f64 - x.ln()*(-x).ln_1p()
        } else if x == 1.0_f64 {
            pi*pi/6.0_f64
        } else if x < 2.0_f64 {
            let l = x.ln();
            li2_approx_f64(1.0_f64 - 1.0_f64/x) + pi*pi/6.0_f64 - l*((1.0_f64 - 1.0_f64/x).ln() + 0.5_f64*l)
        } else {
            let l = x.ln();
            -li2_approx_f64(1.0_f64/x) + pi*pi/3.0_f64 - 0.5_f64*l*l
        }
    }
}

impl Li2<Complex<f64>> for Complex<f64> {
    /// Returns the dilogarithm of a complex number of type
    /// `Complex<f64>`.
    ///
    /// This function has been translated from the
    /// [SPheno](https://spheno.hepforge.org/) package.
    ///
    /// # Example:
    /// ```
    /// use num::complex::Complex;
    /// use polylog::Li2;
    ///
    /// assert!((Complex::new(1.0_f64, 1.0_f64).li2() - Complex::new(0.6168502750680849_f64, 1.4603621167531195_f64)).norm() < std::f64::EPSILON);
    /// ```
    fn li2(&self) -> Complex<f64> {
        let pi = std::f64::consts::PI;

        // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
        // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 19}]
        let bf = [
            - 1.0_f64/4.0_f64,
              1.0_f64/36.0_f64,
            - 1.0_f64/3600.0_f64,
              1.0_f64/211680.0_f64,
            - 1.0_f64/10886400.0_f64,
              1.0_f64/526901760.0_f64,
            - 4.0647616451442255e-11_f64,
              8.9216910204564526e-13_f64,
            - 1.9939295860721076e-14_f64,
              4.5189800296199182e-16_f64,
        ];

        let rz = self.re;
        let iz = self.im;

        if iz == 0.0_f64 {
            if rz <= 1.0_f64 {
                Complex::new(rz.li2(), 0.0_f64)
            } else { // rz > 1
                Complex::new(rz.li2(), -pi*rz.ln())
            }
        } else {
            let nz = self.norm_sqr();

            if nz < std::f64::EPSILON {
                self*(1.0_f64 + 0.25_f64*self)
            } else {
                let (u, rest, sgn) = if rz <= 0.5_f64 {
                    if nz > 1.0_f64 {
                        let l = (-self).cln();
                        (-(1.0_f64 - 1.0_f64/self).cln(), -0.5_f64*l*l - pi*pi/6.0_f64, -1.0_f64)
                    } else { // nz <= 1
                        (-(1.0_f64 - self).cln(), Complex::new(0.0_f64, 0.0_f64), 1.0_f64)
                    }
                } else { // rz > 0.5
                    if nz <= 2.0_f64*rz {
                        let l = -(self).cln();
                        (l, l*(1.0_f64 - self).cln() + pi*pi/6.0_f64, -1.0_f64)
                    } else { // nz > 2*rz
                        let l = (-self).cln();
                        (-(1.0_f64 - 1.0_f64/self).cln(), -0.5_f64*l*l - pi*pi/6.0_f64, -1.0_f64)
                    }
                };

                let u2 = u*u;
                let u4 = u2*u2;
                let sum =
                    u +
                    u2*(bf[0] +
                    u *(bf[1] +
                    u2*(
                        bf[2] +
                        u2*bf[3] +
                        u4*(bf[4] + u2*bf[5]) +
                        u4*u4*(bf[6] + u2*bf[7] + u4*(bf[8] + u2*bf[9]))
                    )));

                sgn*sum + rest
            }
        }
    }
}
