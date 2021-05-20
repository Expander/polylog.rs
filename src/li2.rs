use num::complex::Complex;
use cln::CLn;

/// Provides the dilogarithm function `li2()` of a number of type `T`.
pub trait Li2<T> {
    fn li2(&self) -> T;
}

impl Li2<f64> for f64 {
    /// Returns the real dilogarithm of a real number of type `f64`.
    ///
    /// Implemented as an economized Pade approximation with a maximum
    /// error of 4.16e-18.
    ///
    /// # Example:
    /// ```
    /// use polylog::Li2;
    ///
    /// let z = 1.0;
    /// println!("Li2({}) = {}", z, z.li2());
    /// ```
    fn li2(&self) -> f64 {
        let pi = std::f64::consts::PI;
        let cp = [
            1.0706105563309304277e+0,
           -4.5353562730201404017e+0,
            7.4819657596286408905e+0,
           -6.0516124315132409155e+0,
            2.4733515209909815443e+0,
           -4.6937565143754629578e-1,
            3.1608910440687221695e-2,
           -2.4630612614645039828e-4
        ];
        let cq = [
            1.0000000000000000000e+0,
           -4.5355682121856044935e+0,
            8.1790029773247428573e+0,
           -7.4634190853767468810e+0,
            3.6245392503925290187e+0,
           -8.9936784740041174897e-1,
            9.8554565816757007266e-2,
           -3.2116618742475189569e-3
        ];

        let x = *self;

        // transform to [0, 1/2)
        let (y, r, s) = if x < -1. {
            let l = (1. - x).ln();
            (1./(1. - x), -pi*pi/6. + l*(0.5*l - (-x).ln()), 1.)
        } else if x == -1. {
            return -pi*pi/12.;
        } else if x < 0. {
            let l = (-x).ln_1p();
            (x/(x - 1.), -0.5*l*l, -1.)
        } else if x == 0. {
            return 0.;
        } else if x < 0.5 {
            (x, 0., 1.)
        } else if x < 1. {
            (1. - x, pi*pi/6. - x.ln()*(1. - x).ln(), -1.)
        } else if x == 1. {
            return pi*pi/6.;
        } else if x < 2. {
            let l = x.ln();
            (1. - 1./x, pi*pi/6. - l*((1. - 1./x).ln() + 0.5*l), 1.)
        } else {
            let l = x.ln();
            (1./x, pi*pi/3. - 0.5*l*l, -1.)
        };

        let z = y - 0.25;
        let z2 = z*z;
        let z4 = z2*z2;

        let p = cp[0] + z * cp[1] + z2 * (cp[2] + z * cp[3]) +
                z4 * (cp[4] + z * cp[5] + z2 * (cp[6] + z * cp[7]));
        let q = cq[0] + z * cq[1] + z2 * (cq[2] + z * cq[3]) +
                z4 * (cq[4] + z * cq[5] + z2 * (cq[6] + z * cq[7]));

        r + s*y*p/q
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
    /// extern crate num;
    /// extern crate polylog;
    /// use num::complex::Complex;
    /// use polylog::Li2;
    ///
    /// fn main() {
    ///     let z = Complex::new(1.0, 1.0);
    ///     println!("Li2({}) = {}", z, z.li2());
    /// }
    /// ```
    fn li2(&self) -> Complex<f64> {
        let pi = std::f64::consts::PI;

        // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
        // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 19}]
        let bf = [
            - 1./4.,
              1./36.,
            - 1./3600.,
              1./211680.,
            - 1./10886400.,
              1./526901760.,
            - 4.0647616451442255e-11,
              8.9216910204564526e-13,
            - 1.9939295860721076e-14,
              4.5189800296199182e-16,
        ];

        let rz = self.re;
        let iz = self.im;

        // special cases
        if iz == 0. {
            if rz <= 1. {
                return Complex::new(rz.li2(), 0.0)
            } else { // rz > 1.
                return Complex::new(rz.li2(), -pi*rz.ln())
            }
        }

        let nz = self.norm_sqr();

        if nz < std::f64::EPSILON {
            return self*(1.0 + 0.25*self);
        }

        let (u, rest, sgn) = if rz <= 0.5 {
            if nz > 1. {
                let l = (-self).cln();
                (-(1. - 1. / self).cln(), -0.5 * l * l - pi * pi / 6., -1.)
            } else { // nz <= 1.
                (-(1. - self).cln(), Complex::new(0.,0.), 1.)
            }
        } else { // rz > 0.5
            if nz <= 2.0*rz {
                let l = -(self).cln();
                (l, l * (1. - self).cln() + pi * pi / 6., -1.)
            } else { // nz > 2.0*rz
                let l = (-self).cln();
                (-(1. - 1. / self).cln(), -0.5 * l * l - pi * pi / 6., -1.)
            }
        };

        let u2 = u*u;
        let u4 = u2*u2;
        let sum =
            u +
            u2 * (bf[0] +
            u  * (bf[1] +
            u2 * (
                bf[2] +
                u2*bf[3] +
                u4*(bf[4] + u2*bf[5]) +
                u4*u4*(bf[6] + u2*bf[7] + u4*(bf[8] + u2*bf[9]))
            )));

        sgn * sum + rest
    }
}
