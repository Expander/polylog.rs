use num::complex::Complex;
mod eta;
mod fac;
mod harmonic;
mod zeta;
mod cli;
mod rli;

/// Provides the n-th order polylogarithm function `li()` of a number of type `T`.
pub trait Li<T> {
    fn li(&self, n: i32) -> T;
}

impl Li<Complex<f64>> for Complex<f64> {
    /// Returns the complex n-th order polylogarithm of a complex
    /// number of type `Complex<f64>` for all integers `n`.
    ///
    /// The implementation for `n < 0` is an adaption of
    /// [[arxiv:2010.09860]].
    ///
    /// [arxiv:2010.09860]: https://arxiv.org/abs/2010.09860
    ///
    /// # Example:
    /// ```
    /// use num::complex::Complex;
    /// use polylog::Li;
    ///
    /// let z = Complex::new(1.0, 1.0);
    /// let n = 10;
    /// println!("Li({},{}) = {}", n, z, z.li(n));
    /// ```
    fn li(&self, n: i32) -> Complex<f64> {
        cli::cli(n, *self)
    }
}

impl Li<f64> for f64 {
    /// Returns the real n-th order polylogarithm of a real number of
    /// type `f64` for all integers `n`.
    ///
    /// The implementation for `n < 0` is an adaption of
    /// [[arxiv:2010.09860]].
    ///
    /// [arxiv:2010.09860]: https://arxiv.org/abs/2010.09860
    ///
    /// # Example:
    /// ```
    /// use polylog::Li;
    ///
    /// let z = 1.0;
    /// let n = 10;
    /// println!("Li({},{}) = {}", n, z, z.li(n));
    /// ```
    fn li(&self, n: i32) -> f64 {
        rli::rli(n, *self)
    }
}
