use num::complex::Complex;
use crate::cln::CLn;

/// Provides the 4-th order polylogarithm function `li4()` of a
/// number of type `T`.
pub trait Li4<T> {
    fn li4(&self) -> T;
}

impl Li4<f64> for f64 {
    /// Returns the fourth order polylogarithm of a real number of type `f64`.
    ///
    /// # Example:
    /// ```
    /// use polylog::Li4;
    ///
    /// assert!((1.0_f64.li4() - 1.0823232337111382_f64).abs() < std::f64::EPSILON);
    /// ```
    fn li4(&self) -> f64 {
        let z2 = 1.6449340668482264;
        let z4 = 1.0823232337111382;
        let x = *self;

        // transform x to y in [-1,1]
        let (y, rest, sgn) = if x < -1.0 {
            let l = (-x).ln();
            let l2 = l*l;
            (1.0/x, -7.0/4.0*z4 + l2*(-0.5*z2 - 1.0/24.0*l2), -1.0)
        } else if x == -1.0 {
            return -7.0/8.0*z4
        } else if x < 1.0 {
            (x, 0.0, 1.0)
        } else if x == 1.0 {
            return z4
        } else { // x > 1.0
            let l = x.ln();
            let l2 = l*l;
            (1.0/x, 2.0*z4 + l2*(z2 - 1.0/24.0*l2), -1.0)
        };

        if y < 0.0 {
            sgn*li4_neg(y) + rest
        } else if y < 0.5 {
            sgn*li4_half(y) + rest
        } else if y < 0.8 {
            sgn*li4_mid(y) + rest
        } else { // y <= 1.0
            sgn*li4_one(y) + rest
        }
    }
}

// Li_4(x) for x in [-1,0]
fn li4_neg(x: f64) -> f64 {
    let cp = [
        0.9999999999999999952e+0, -1.8532099956062184217e+0,
        1.1937642574034898249e+0, -3.1817912243893560382e-1,
        3.2268284189261624841e-2, -8.3773570305913850724e-4
    ];
    let cq = [
        1.0000000000000000000e+0, -1.9157099956062165688e+0,
        1.3011504531166486419e+0, -3.7975653506939627186e-1,
        4.5822723996558783670e-2, -1.8023912938765272341e-3,
        1.0199621542882314929e-5
    ];

    let x2 = x*x;
    let x4 = x2*x2;
    let p = cp[0] + x * cp[1] + x2 * (cp[2] + x * cp[3]) +
            x4 * (cp[4] + x * cp[5]);
    let q = cq[0] + x * cq[1] + x2 * (cq[2] + x * cq[3]) +
            x4 * (cq[4] + x * cq[5] + x2 * cq[6]);

    x*p/q
}

// Li_4(x) for x in [0,1/2]
fn li4_half(x: f64) -> f64 {
    let cp = [
        1.0000000000000000414e+0, -2.0588072418045364525e+0,
        1.4713328756794826579e+0, -4.2608608613069811474e-1,
        4.2975084278851543150e-2, -6.8314031819918920802e-4
    ];
    let cq = [
        1.0000000000000000000e+0, -2.1213072418045207223e+0,
        1.5915688992789175941e+0, -5.0327641401677265813e-1,
        6.1467217495127095177e-2, -1.9061294280193280330e-3
    ];

    let x2 = x*x;
    let x4 = x2*x2;
    let p = cp[0] + x * cp[1] + x2 * (cp[2] + x * cp[3]) +
            x4 * (cp[4] + x * cp[5]);
    let q = cq[0] + x * cq[1] + x2 * (cq[2] + x * cq[3]) +
            x4 * (cq[4] + x * cq[5]);

    x*p/q
}

// Li_4(x) for x in [1/2,8/10]
fn li4_mid(x: f64) -> f64 {
    let cp = [
        3.2009826406098890447e-9, 9.9999994634837574160e-1,
       -2.9144851228299341318e+0, 3.1891031447462342009e+0,
       -1.6009125158511117090e+0, 3.5397747039432351193e-1,
       -2.5230024124741454735e-2
    ];
    let cq = [
        1.0000000000000000000e+0, -2.9769855248411488460e+0,
        3.3628208295110572579e+0, -1.7782471949702788393e+0,
        4.3364007973198649921e-1, -3.9535592340362510549e-2,
        5.7373431535336755591e-4
    ];

    let x2 = x*x;
    let x4 = x2*x2;
    let p = cp[0] + x * cp[1] + x2 * (cp[2] + x * cp[3]) +
            x4 * (cp[4] + x * cp[5] + x2 * cp[6]);
    let q = cq[0] + x * cq[1] + x2 * (cq[2] + x * cq[3]) +
            x4 * (cq[4] + x * cq[5] + x2 * cq[6]);

    p/q
}

// Li_4(x) for x in [8/10,1]
fn li4_one(x: f64) -> f64 {
    let z2 = 1.6449340668482264;
    let z3 = 1.2020569031595943;
    let z4 = 1.0823232337111382;
    let l = x.ln();
    let l2 = l*l;

    z4 +
    l*(z3 +
    l*(0.5*z2 +
    l*(11.0/36.0 - 1.0/6.0*(-l).ln() +
    l*(-1.0/48.0 +
    l*(-1.0/1440.0 +
    l2*(1.0/604800.0 - 1.0/91445760.0*l2))))))
}

impl Li4<Complex<f64>> for Complex<f64> {
    /// Returns the fourth order polylogarithm of a complex number of type
    /// `Complex<f64>`.
    ///
    /// # Example:
    /// ```
    /// use num::complex::Complex;
    /// use polylog::Li4;
    ///
    /// assert!((Complex::new(1.0_f64, 1.0_f64).li4() - Complex::new(0.9593189135784193_f64, 1.1380391966769828_f64)).norm() < 2.0_f64*std::f64::EPSILON);
    /// ```
    fn li4(&self) -> Complex<f64> {
        let pi  = std::f64::consts::PI;
        let pi2 = pi*pi;
        let z4  = 1.0823232337111382;

        if self.im == 0.0 {
            if self.re <= 1.0 {
                return Complex::new(self.re.li4(), 0.0)
            } else { // rz > 1.0
                let l = self.re.ln();
                return Complex::new(self.re.li4(), -pi/6.0*l*l*l)
            }
        }

        let nz  = self.norm();
        let pz  = self.arg();
        let lnz = nz.ln();

        if lnz*lnz + pz*pz < 1.0 { // |log(z)| < 1
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

        if nz <= 1.0 {
            cli4_unit_circle(-(1.0 - self).cln())
        } else { // nz > 1.0
            let pi4  = pi2*pi2;
            let arg = if pz > 0.0 { pz - pi } else { pz + pi };
            let lmz = Complex::new(lnz, arg); // (-self).cln()
            let lmz2 = lmz*lmz;
            -cli4_unit_circle(-(1.0 - 1.0/self).cln()) + 1.0/360.0*(-7.0*pi4 + lmz2*(-30.0*pi2 - 15.0*lmz2))
        }
    }
}

/// series approximation of Li3(z) for |z| <= 1
/// in terms of x = -ln(1 - z)
fn cli4_unit_circle(x: Complex<f64>) -> Complex<f64> {
    let bf  = [
        1.0                   , -7.0/16.0              ,
        1.1651234567901235e-01, -1.9820601851851852e-02,
        1.9279320987654321e-03, -3.1057098765432099e-05,
       -1.5624009114857835e-05,  8.4851235467732066e-07,
        2.2909616603189711e-07, -2.1832614218526917e-08,
       -3.8828248791720156e-09,  5.4462921032203321e-10,
        6.9608052106827254e-11, -1.3375737686445215e-11,
       -1.2784852685266572e-12,  3.2605628580248922e-13,
        2.3647571168618257e-14, -7.9231351220311617e-15,
    ];

    let x2 = x*x;
    let x4 = x2*x2;
    let x8 = x4*x4;

    x*bf[0] +
    x2*(bf[1] + x*bf[2]) +
    x4*(bf[3] + x*bf[4] + x2*(bf[5] + x*bf[6])) +
    x8*(bf[7] + x*bf[8] + x2*(bf[9] + x*bf[10]) +
        x4*(bf[11] + x*bf[12] + x2*(bf[13] + x*bf[14]))) +
    x8*x8*(bf[15] + x*bf[16] + x2*bf[17])
}
