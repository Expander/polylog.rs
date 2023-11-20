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
    /// let x = 2.0;
    /// println!("Li0({}) = {}", x, x.li0());
    /// ```
    fn li0(&self) -> f64 {
        self/(1.0 - self)
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
    /// let z = Complex::new(1.0, 1.0);
    /// println!("Li0({}) = {}", z, z.li0());
    /// ```
    fn li0(&self) -> Complex<f64> {
        self/(1.0 - self)
    }
}
