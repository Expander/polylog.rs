use num::complex::Complex;

/// Provides an implementation of the complex logarithm `cln()` of a
/// number of type `T`, where the imaginary part of the logarithm is
/// within (-pi, pi]. In particular, arguments z with Im(z) == -0.0
/// are treated as having Im(z) == 0.0.
///
/// Note, that this behavior is different from the complex ln(), where
/// the imaginary part is within [-pi, pi].
pub trait CLn<T> {
    fn cln(&self) -> T;
}

impl CLn<Complex<f64>> for Complex<f64> {
    fn cln(&self) -> Complex<f64> {
        let z = Complex::new(
            self.re,
            // convert -0.0 to 0.0
            if self.im == 0.0 { 0.0 } else { self.im },
        );
        Complex::new(z.norm().ln(), z.arg())
    }
}

#[test]
fn test_cln() {
    let x = Complex::new(1.0, 0.0);

    assert!((-x).cln() == Complex::new(0.0, std::f64::consts::PI));
    assert!((-x).ln() == Complex::new(0.0, -std::f64::consts::PI));
}
