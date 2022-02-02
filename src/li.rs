use {Li0, Li1, Li2, Li3, Li4};

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
    /// println!("Li({},{}) = {}", n, z, z.li(n));
    /// ```
    fn li(&self, n: i64) -> f64 {
        if n < 0 {
            panic!("li(n) not implemented for n < 0 (given value: n = {})", n);
        } else if n == 0 {
            return self.li0();
        } else if n == 1 {
            return self.li1();
        } else if n == 2 {
            return self.li2();
        } else if n == 3 {
            return self.li3();
        } else if n == 4 {
            return self.li4();
        }

        0.0
    }
}
