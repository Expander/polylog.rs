use num::complex::Complex;
use num::{Float, FromPrimitive};

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

impl<T: Float + FromPrimitive> CLn<Complex<T>> for Complex<T> {
    fn cln(&self) -> Complex<T> {
        if self.im == T::zero() && self.re > T::zero() {
            Complex::new(self.re.ln(), T::zero())
        } else if self.im == T::zero() {
            Complex::new((-self.re).ln(), T::from_f64(3.1415926535897932_f64).unwrap())
        } else {
            self.ln()
        }
    }
}

#[test]
fn test_cln() {
    let x32 = Complex::new(1.0_f32, 0.0_f32);
    let x64 = Complex::new(1.0_f64, 0.0_f64);

    assert!(((-x32).cln() - Complex::new(0.0_f32,  std::f32::consts::PI)).norm() < 4.0_f32*std::f32::EPSILON);
    assert!(((-x32).ln()  - Complex::new(0.0_f32, -std::f32::consts::PI)).norm() < 4.0_f32*std::f32::EPSILON);
    assert!((-x64).cln() == Complex::new(0.0_f64,  std::f64::consts::PI));
    assert!((-x64).ln()  == Complex::new(0.0_f64, -std::f64::consts::PI));
}
