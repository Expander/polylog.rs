use num::complex::Complex;
use cln::CLn;

/// Provides the fourth order polylogarithm function `li4()` of a
/// number of type `T`.
pub trait Li4<T> {
    fn li4(&self) -> T;
}

impl Li4<Complex<f64>> for Complex<f64> {
    /// Returns the fourth order polylogarithm of a complex number of type
    /// `Complex<f64>`.
    ///
    /// # Example:
    /// ```
    /// extern crate num;
    /// extern crate polylog;
    /// use num::complex::Complex;
    /// use polylog::Li4;
    ///
    /// fn main() {
    ///     let z = Complex::new(1.0, 1.0);
    ///     println!("Li4({}) = {}", z, z.li4());
    /// }
    /// ```
    fn li4(&self) -> Complex<f64> {
        let pi  = std::f64::consts::PI;
        let pi2 = pi*pi;
        let z4  = 1.0823232337111382;
        let bf  = [
            1.                    , -7./16.                ,
            1.1651234567901235e-01, -1.9820601851851852e-02,
            1.9279320987654321e-03, -3.1057098765432099e-05,
           -1.5624009114857835e-05,  8.4851235467732066e-07,
            2.2909616603189711e-07, -2.1832614218526917e-08,
           -3.8828248791720156e-09,  5.4462921032203321e-10,
            6.9608052106827254e-11, -1.3375737686445215e-11,
           -1.2784852685266572e-12,  3.2605628580248922e-13,
            2.3647571168618257e-14, -7.9231351220311617e-15,
        ];

        if self.im == 0.0 {
            if self.re == 0.0 {
                return Complex::new(0., 0.);
            }
            if self.re == 1.0 {
                return Complex::new(z4, 0.);
            }
            if self.re == -1.0 {
                return Complex::new(-7./8.*z4, 0.);
            }
        }

        let nz  = self.norm_sqr();
        let pz  = self.arg();
        let lnz = 0.5*nz.ln();

        if lnz*lnz + pz*pz < 1. { // |log(z)| < 1
            let u  = Complex::new(lnz, pz);
            let u2 = u*u;
            let u4 = u2*u2;
            let u8 = u4*u4;
            let c1 = 1.2020569031595943; // zeta(3)
            let c2 = 0.82246703342411322;
            let c3 = (11.0/6.0 - (-u).cln())/6.0;
            let c4 = -1.0/48.0;

            let cs = [
                -6.9444444444444444e-04, 1.6534391534391534e-06,
                -1.0935444136502338e-08, 1.0438378493934049e-10,
                -1.2165942300622435e-12, 1.6130006528350101e-14,
                -2.3428810452879340e-16
            ];

            return z4 + u2 * (c2 + u2 * c4) +
                u * (
                    c1 +
                    c3*u2 +
                    u4*(cs[0] + u2*cs[1]) +
                    u8*(cs[2] + u2*cs[3] + u4*(cs[4] + u2*cs[5])) +
                    u8*u8*cs[6]
                );
        }

        let (u, rest, sgn) = if nz <= 1.0 {
            (-(1. - self).cln(), Complex::new(0.,0.), 1.)
        } else { // nz > 1.0
            let pi4  = pi2*pi2;
            let arg = if pz > 0.0 { pz - pi } else { pz + pi };
            let lmz = Complex::new(lnz, arg); // (-self).cln()
            let lmz2 = lmz*lmz;
            (-(1. - 1./self).cln(), 1./360.*(-7.*pi4 + lmz2*(-30.*pi2 - 15.*lmz2)), -1.)
        };

        let u2 = u*u;
        let u4 = u2*u2;
        let u8 = u4*u4;

        rest + sgn * (
           u*bf[0] +
           u2*(bf[1] + u*bf[2]) +
           u4*(bf[3] + u*bf[4] + u2*(bf[5] + u*bf[6])) +
           u8*(bf[7] + u*bf[8] + u2*(bf[9] + u*bf[10]) +
               u4*(bf[11] + u*bf[12] + u2*(bf[13] + u*bf[14]))) +
           u8*u8*(bf[15] + u*bf[16] + u2*bf[17])
        )
    }
}
