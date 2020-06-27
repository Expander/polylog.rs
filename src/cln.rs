use num::complex::Complex;

/// Provides an implementation of the complex logarithm `cln()` of a
/// number of type `T`, where the imaginary part of the logarithm is
/// within (-Pi,Pi].
pub trait CLn<T> {
    fn cln(&self) -> T;
}

impl CLn<Complex<f64>> for Complex<f64> {
    fn cln(&self) -> Complex<f64> {
        let z = Complex::new(
            self.re,
            // convert -0.0 to 0.0
            if self.im == 0. { 0. } else { self.im },
        );
        Complex::new(0.5*z.norm_sqr().ln(), z.arg())
    }
}
