use std;
use num::complex::Complex;
use num::Float;

/// Provides the trilogarithm function `li3()` of a number of type
/// `T`.
pub trait Li3<T> {
    fn li3(&self) -> T;
}

/// Returns the trilogarithm of a real number of type `f64`.
///
/// # Example:
/// ```
/// extern crate num;
/// extern crate polylog;
/// use num::complex::Complex;
/// use polylog::Li3;
///
/// fn main() {
///     let z = Complex::new(1.0, 1.0);
///     println!("Li3({}) = {}", z, z.li3());
/// }
/// ```
impl Li3<Complex<f64>> for Complex<f64> {
    fn li3(&self) -> Complex<f64> {
        let pi  = 3.1415926535897932384626433832795;
        let pi2 = pi*pi;
        let eps = std::f64::EPSILON;
        let z3  = 1.2020569031595942853997381615114;
        let bf  = vec![
            1., -3./8., 17./216., -5./576.,
             0.00012962962962962962962962962962963,  0.000081018518518518518518518518518519,
            -3.4193571608537594932152755282007e-06, -1.3286564625850340136054421768707e-06 ,
             8.6608717561098513479465860418241e-08,  2.5260875955320399764844209288654e-08 ,
            -2.1446944683640647609338850757365e-09, -5.1401106220129789153358176927200e-10 ,
             5.2495821146008294363940888085581e-11,  1.0887754406636318375372971570425e-11 ,
            -1.2779396094493695305581831754072e-12, -2.3698241773087452099797778810124e-13 ,
             3.1043578879654622942847532704656e-14,  5.2617586299125060841318392511225e-15 ,
            -7.5384795499492653659925014322677e-16, -1.1862322577752285253082500951246e-16 ,
             1.8316979965491383382089273121282e-17,  2.7068171031837350151490734712617e-18 ,
            -4.4554338978296388264326309921763e-19, -6.2375484922556946503653222473984e-20 ,
             1.0851521534874534913136560996864e-20,  1.4491174866036081930734904966528e-21 ,
            -2.6466339754458990334740891186144e-22, -3.3897653488510104721925816586081e-23 ,
             6.4640477336033108890325309821953e-24,  7.9758344896024124242092227259050e-25 ,
            -1.5809178790287483355921117629383e-25, -1.8861499729622868193110225398853e-26 ,
             3.8715536638418473303997127188831e-27,  4.4801175002345607304865389832051e-28 ,
            -9.4930338719118361264175367602770e-29, -1.0682813809077381224018214303381e-29 ,
             2.3304478936103051860078519901928e-30,  2.5560775726519754080563569828670e-31 ,
            -5.7274216061372596844727445803306e-32, -6.1347132137964235825854929689777e-33
        ];

        if is_close(self, 0., eps) {
            return Complex::new(0., 0.);
        }
        if is_close(self, 1., eps) {
            return Complex::new(z3, 0.);
        }
        if is_close(self, 1., 0.02) {
            let i    = Complex::i();
            let ipi  = Complex::new(0., pi);
            let zm1  = *self - 1.;
            let lzm1 = zm1.cln();
            let ceil = if zm1.arg() > 0. { 1. } else { 0. };

            return z3
                + pi2*zm1/6.
                + (ceil*ipi - lzm1/2. - 1./12.*pow2(3.*i + pi))*pow2(zm1)
                + (lzm1/2. + 1./36.*(-21. - 18.*(-1. + 2.*ceil)*ipi + 2.*pi2))*pow3(zm1)
                + (-11./24.*lzm1 + 1./288.*(131. + 132.*(-1. + 2.*ceil)*ipi - 12.*pi2))*pow4(zm1)
                + (5./12.*lzm1 + 1./720.*(-265. - 300.*(-1. + 2.*ceil)*ipi + 24.*pi2))*pow5(zm1)
                + (-137./360.*lzm1 + 1./7200.*(2213. + 2740.*(-1. + 2.*ceil)*ipi - 200.*pi2))*pow6(zm1)
                + (-947./3600. - 7./20.*(-1. + 2.*ceil)*ipi + 7./20.*lzm1 + pi2/42.)*pow7(zm1)
                + (647707./2822400. + 363.*(-1. + 2.*ceil)*ipi/1120. - 363./1120.*lzm1 - pi2/48.)*pow8(zm1);
        }
        if is_close(self, -1., eps) {
            return Complex::new(-0.75*z3, 0.);
        }
        if is_close(self, 0.5, eps) {
            let ln2  = (2.).ln();
            let ln23 = ln2.powi(3);
            return Complex::new((-2.*pi2*ln2 + 4.*ln23 + 21.*z3)/24., 0.);
        }

        let (u, mut sum) = if self.norm() <= 1. {
            (-(1. - self).cln(), Complex::new(0.,0.))
        } else { // az > 1.
            (-(1. - 1./self).cln(), -pow3((-self).cln())/6. - pi2/6.*(-self).cln())
        };

        let mut p = Complex::new(1.,0.);

        for b in bf {
            p *= u;
            sum += b*p;
        }

        sum
    }
}

fn is_close(a : &Complex<f64>, b : f64, eps : f64) -> bool {
    (a.re - b).abs() < eps && (a.im).abs() < eps
}

fn pow2(z : Complex<f64>) -> Complex<f64> {
    z * z
}

fn pow3(z : Complex<f64>) -> Complex<f64> {
    z * z * z
}

fn pow4(z : Complex<f64>) -> Complex<f64> {
    z * z * z * z
}

fn pow5(z : Complex<f64>) -> Complex<f64> {
    z * z * z * z * z
}

fn pow6(z : Complex<f64>) -> Complex<f64> {
    z * z * z * z * z * z
}

fn pow7(z : Complex<f64>) -> Complex<f64> {
    z * z * z * z * z * z * z
}

fn pow8(z : Complex<f64>) -> Complex<f64> {
    z * z * z * z * z * z * z * z
}

trait CLn<T> {
    fn cln(&self) -> T;
}

impl CLn<Complex<f64>> for Complex<f64> {
    fn cln(&self) -> Complex<f64> {
        Complex::new(
            if self.re == 0. { 0. } else { self.re },
            if self.im == 0. { 0. } else { self.im },
        ).ln()
    }
}
