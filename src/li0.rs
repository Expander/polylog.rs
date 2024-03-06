use num::complex::Complex;

/// Provides the 0-th order polylogarithm function `li0()` of a
/// number of type `T`.
pub trait Li0<T> {
    fn li0(&self) -> T;
}

impl Li0<f64> for f64 {
    /// Returns the real 0th order polylogarithm of a real number of
    /// type `f64`.
    ///
    /// # Example:
    /// ```
    /// use polylog::Li0;
    ///
    /// assert!((2.0_f64.li0() + 2.0_f64).abs() < std::f64::EPSILON);
    /// ```
    fn li0(&self) -> f64 {
        if *self == 0.0 {
            *self
        } else {
            self/(1.0 - self)
        }
    }
}

impl Li0<Complex<f64>> for Complex<f64> {
    /// Returns the 0th order polylogarithm of a complex number of
    /// type `Complex<f64>`.
    ///
    /// # Example:
    /// ```
    /// use num::complex::Complex;
    /// use polylog::Li0;
    ///
    /// assert!((Complex::new(1.0_f64, 1.0_f64).li0() - Complex::new(-1.0_f64, 1.0_f64)).norm() < std::f64::EPSILON);
    /// ```
    fn li0(&self) -> Complex<f64> {
        if *self == Complex::new(0.0, 0.0) {
            *self
        } else {
            self/(1.0 - self)
        }
    }
}
