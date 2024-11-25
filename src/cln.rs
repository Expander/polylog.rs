use num::complex::Complex;
use num::Float;

trait Pi<T> {
    fn pi() -> T;
}

impl Pi<f32> for f32 {
    fn pi() -> f32 { std::f32::consts::PI }
}

impl Pi<f64> for f64 {
    fn pi() -> f64 { std::f64::consts::PI }
}

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

impl<T: Float + Pi<T>> CLn<Complex<T>> for Complex<T> {
    fn cln(&self) -> Complex<T> {
        if self.im == T::zero() && self.re > T::zero() {
            Complex::new(self.re.ln(), T::zero())
        } else if self.im == T::zero() && self.re < T::zero() {
            Complex::new((-self.re).ln(), T::pi())
        } else {
            self.ln()
        }
    }
}

#[test]
fn test_cln() {
    // test positive zero
    let cnegone32 = Complex::new(-1.0_f32, 0.0_f32);
    let cnegone64 = Complex::new(-1.0_f64, 0.0_f64);

    assert!((cnegone32.cln() - Complex::new(0.0_f32, std::f32::consts::PI)).norm() < 4.0_f32*std::f32::EPSILON);
    assert!((cnegone32.ln()  - Complex::new(0.0_f32, std::f32::consts::PI)).norm() < 4.0_f32*std::f32::EPSILON);
    assert!(cnegone64.cln() == Complex::new(0.0_f64, std::f64::consts::PI));
    assert!(cnegone64.ln()  == Complex::new(0.0_f64, std::f64::consts::PI));

    // test negative zero
    let cone32 = Complex::new(1.0_f32, 0.0_f32);
    let cone64 = Complex::new(1.0_f64, 0.0_f64);

    assert!(((-cone32).cln() - Complex::new(0.0_f32,  std::f32::consts::PI)).norm() < 4.0_f32*std::f32::EPSILON);
    assert!(((-cone32).ln()  - Complex::new(0.0_f32, -std::f32::consts::PI)).norm() < 4.0_f32*std::f32::EPSILON);
    assert!((-cone64).cln() == Complex::new(0.0_f64,  std::f64::consts::PI));
    assert!((-cone64).ln()  == Complex::new(0.0_f64, -std::f64::consts::PI));

    // test zero input
    let rzero32 = 0.0_f32;
    let rzero64 = 0.0_f64;

    assert!(Complex::new( rzero32, rzero32).cln() == Complex::new(std::f32::NEG_INFINITY, 0.0));
    assert!(Complex::new(-rzero32, rzero32).cln().re == std::f32::NEG_INFINITY);
    assert!((Complex::new(-rzero32, rzero32).cln().im - std::f32::consts::PI).abs() < 4.0_f32*std::f32::EPSILON);
    assert!(Complex::new( rzero64, rzero64).cln() == Complex::new(std::f64::NEG_INFINITY, 0.0));
    assert!(Complex::new(-rzero64, rzero64).cln() == Complex::new(std::f64::NEG_INFINITY, std::f64::consts::PI));
}
