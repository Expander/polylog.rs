use num::complex::Complex;
use crate::cln::CLn;

/// Provides the 1st order polylogarithm function `li1()` of a
/// number of type `T`.
pub trait Li1<T> {
    fn li1(&self) -> T;
}

impl Li1<f64> for f64 {
    /// Returns the real first order polylogarithm of a real number of
    /// type `f64`.
    ///
    /// # Example:
    /// ```
    /// use polylog::Li1;
    ///
    /// let x = 2.0;
    /// println!("Li1({}) = {}", x, x.li1());
    /// ```
    fn li1(&self) -> f64 {
        let x = *self;
        if x < 1.0 {
            -(-x).ln_1p()
        } else if x == 1.0 {
            std::f64::INFINITY
        } else { // x > 1.0
            -(x - 1.0).ln()
        }
    }
}

impl Li1<Complex<f64>> for Complex<f64> {
    /// Returns the first order polylogarithm of a complex number of
    /// type `Complex<f64>`.
    ///
    /// # Example:
    /// ```
    /// extern crate num;
    /// use num::complex::Complex;
    /// use polylog::Li1;
    ///
    /// fn main() {
    ///     let z = Complex::new(1.0, 1.0);
    ///     println!("Li1({}) = {}", z, z.li1());
    /// }
    /// ```
    fn li1(&self) -> Complex<f64> {
        -(1.0 - self).cln()
    }
}
