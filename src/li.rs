/// Provides the n-th order polylogarithm function `li()` of a number of type `T`.
pub trait Li<T> {
    fn li(&self, i64) -> T;
}

impl Li<f64> for f64 {
    /// Returns the real n-th order polylogarithm of a real number of
    /// type `f64`.
    ///
    /// # Example:
    /// ```
    /// use polylog::Li;
    ///
    /// let z = 1.0;
    /// let n = 10;
    /// println!("Li({},{}) = {}", n, z, z.li());
    /// ```
    fn li(&self, n: i64) -> f64 {
        0.0
    }
}
