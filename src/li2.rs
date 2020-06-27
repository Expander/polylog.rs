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
        let p = horner(z, &cp);
        let q = horner(z, &cq);

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
        let nz = self.norm_sqr();

        // special cases
        if iz == 0. {
            if rz <= 1. {
                return Complex::new(rz.li2(), 0.0)
            } else { // rz > 1.
                return Complex::new(rz.li2(), -pi*rz.ln())
            }
        } else if nz < std::f64::EPSILON {
            return *self;
        }

        let (cy, cz, sgn) = if rz <= 0.5 {
            if nz > 1. {
                let l = (-self).cln();
                (-0.5 * l * l - pi * pi / 6., -(1. - 1. / self).cln(), -1.)
            } else { // nz <= 1.
                (Complex::new(0.,0.), -(1. - self).cln(), 1.)
            }
        } else { // rz > 0.5
            if nz <= 2.0*rz {
                let l = -(self).cln();
                (l * (1. - self).cln() + pi * pi / 6., l, -1.)
            } else { // nz > 2.0*rz
                let l = (-self).cln();
                (-0.5 * l * l - pi * pi / 6., -(1. - 1. / self).cln(), -1.)
            }
        };

        let cz2 = cz*cz;
        let cz4 = cz2*cz2;
        let sum =
            cz +
            cz2 * (bf[0] +
            cz  * (bf[1] +
            cz2 * (
                bf[2] +
                cz2*bf[3] +
                cz4*(bf[4] + cz2*bf[5]) +
                cz4*cz4*(bf[6] + cz2*bf[7] + cz4*(bf[8] + cz2*bf[9]))
            )));

        sgn * sum + cy
    }
}

/// evaluation of polynomial P(x) with coefficients `coeffs`
fn horner(x: f64, coeffs: &[f64]) -> f64 {
    let mut p = 0.;
    for c in coeffs.iter().rev() {
        p = p*x + c;
    }
    p
}
