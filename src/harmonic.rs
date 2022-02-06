/// returns n-th harmonic number, n > 0
pub fn harmonic(n: i32) -> f64 {
    if n <= 0 {
        panic!("harmonic not implemented for n <= 0 (given value: n = {})",  n)
    } else if n < 20 {
        let mut sum = 1.0;
        for k in 2..=n {
            sum += (k as f64).recip();
        }
        sum
    } else {
        let eulergamma = 0.57721566490153286;
        eulergamma + digamma(n + 1)
    }
}

#[test]
fn test_harmonic() {
    let abs = |x: f64| if x >= 0.0 { x } else { -x };

    assert!(harmonic(1) == 1.0);
    assert!(harmonic(2) == 3.0/2.0);
    assert!(harmonic(3) == 11.0/6.0);
    assert!(abs(harmonic( 4) - 25.0/12.0         ) < 1e-14);
    assert!(abs(harmonic(19) - 3.5477396571436819) < 1e-14);
    assert!(abs(harmonic(20) - 3.5977396571436819) < 1e-14);
    assert!(abs(harmonic(21) - 3.6453587047627295) < 1e-14);
}

/// digamma for integer n > 0, following
/// [K.S. Kölbig: Programs for computing the logarithm of the gamma
/// function, and the digamma function, for complex argument, Computer
/// Physics Communications, Volume 4, Issue 2, 1972, Pages 221-226, ISSN
/// 0010-4655, https://doi.org/10.1016/0010-4655(72)90012-4]
fn digamma(n: i32) -> f64 {
    // Table[BernoulliB[2n]/(2 n), {n,1,8}]
    let c = [
        0.083333333333333333, -0.0083333333333333333,  0.0039682539682539683,
       -0.0041666666666666667, 0.0075757575757575758, -0.021092796092796093,
        0.083333333333333333, -0.44325980392156863
    ];

    if n <= 0 {
        panic!("digamma not implemented for n <= 0 (given value: n = {})", n);
    }

    // map potentially small n to m >= 7
    let (m, res) = if n < 7 { // recurrence formula
        let mut res = 0.0;
        for nu in 1..(7 - n) {
            res -= ((n + nu) as f64).recip();
        }
        (7.0, res - (n as f64).recip())
    } else {
        (n as f64, 0.0)
    };

    let t = m.recip();
    let t2 = t*t;

    res + m.ln() - 0.5*t
    - t2*(c[0] + t2*(c[1] + t2*(c[2] + t2*(c[3] + t2*(c[4] + t2*(c[5] + t2*(c[6] + t2*c[7])))))))
}

#[test]
fn test_digamma() {
    let abs = |x: f64| if x >= 0.0 { x } else { -x };

    assert!(abs(digamma( 1) - -0.57721566490153286) < 1e-14);
    assert!(abs(digamma( 2) -  0.42278433509846714) < 1e-14);
    assert!(abs(digamma( 3) -  0.92278433509846714) < 1e-14);
    assert!(abs(digamma(10) -  2.2517525890667211 ) < 1e-15);
}
