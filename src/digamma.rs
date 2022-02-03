/// digamma for integer n > 0, following
/// [K.S. Kölbig: Programs for computing the logarithm of the gamma
/// function, and the digamma function, for complex argument, Computer
/// Physics Communications, Volume 4, Issue 2, 1972, Pages 221-226, ISSN
/// 0010-4655, https://doi.org/10.1016/0010-4655(72)90012-4]
pub fn digamma(n: i32) -> f64 {
    // Table[BernoulliB[2n]/(2 n), {n,1,8}]
    let c = [
        0.083333333333333333, -0.0083333333333333333,  0.0039682539682539683,
       -0.0041666666666666667, 0.0075757575757575758, -0.021092796092796093,
        0.083333333333333333, -0.44325980392156863
    ];

    if n <= 0 {
        panic!("digamma not implemented for n <= 0 (given value: n = {})", n);
    }

    let mut res = 0.0;

    let m = if n < 7 { // recurrence formula
        for nu in 1..(7 - n) {
            res -= 1.0/((n + nu) as f64);
        }
        res -= 1.0/(n as f64);
        7
    } else {
        n
    };

    let t = 1.0/(m as f64);
    let t2 = t*t;

    res + (m as f64).ln() - 0.5*t
    - t2*(c[0] + t2*(c[1] + t2*(c[2] + t2*(c[3] + t2*(c[4] + t2*(c[5] + t2*(c[6] + t2*c[7])))))))
}

#[test]
fn test_values() {
    let abs = |x: f64| if x >= 0.0 { x } else { -x };

    assert!(abs(digamma( 1) - -0.57721566490153286) < 1e-14);
    assert!(abs(digamma( 2) -  0.42278433509846714) < 1e-14);
    assert!(abs(digamma( 3) -  0.92278433509846714) < 1e-14);
    assert!(abs(digamma(10) -  2.2517525890667211 ) < 1e-15);
}
