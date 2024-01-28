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
    /// assert!((2.0_f64.li1()).abs() < std::f64::EPSILON);
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
    /// use num::complex::Complex;
    /// use polylog::Li1;
    ///
    /// assert!((Complex::new(1.0_f64, 1.0_f64).li1() - Complex::new(0.0_f64, 1.5707963267948966_f64)).norm() < std::f64::EPSILON);
    /// ```
    fn li1(&self) -> Complex<f64> {
        -(1.0 - self).cln()
    }
}
