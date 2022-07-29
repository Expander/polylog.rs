Polylog
=======

[![Build Status](https://github.com/Expander/polylog/workflows/test/badge.svg)](https://github.com/Expander/polylog/actions)
[![Documentation](https://docs.rs/polylog/badge.svg)](https://docs.rs/polylog/)

The Polylog package provides Rust implementations of real and complex
polylogarithms, including the dilogarithm and trilogarithm.

The Polylog package depends on the `num` crate.


Example
-------

```rust
extern crate num;
extern crate polylog;
use num::complex::Complex;
use polylog::{Li, Li0, Li1, Li2, Li3, Li4, Li5, Li6};

fn main() {
    let x = 1.0;
    let z = Complex::new(1.0, 1.0);
    let n = 10;

    // real polylogarithms for real arguments
    println!("Li0({}) = {}", x, x.li0());      // Re[Li_0(x)]
    println!("Li1({}) = {}", x, x.li1());      // Re[Li_1(x)]
    println!("Li2({}) = {}", x, x.li2());      // Re[Li_2(x)] (dilogarithm)
    println!("Li3({}) = {}", x, x.li3());      // Re[Li_3(x)] (trilogarithm)
    println!("Li4({}) = {}", x, x.li4());      // Re[Li_4(x)]
    println!("Li_{}({}) = {}", n, x, x.li(n)); // Re[Li_n(x)]

    // complex polylogarithms for complex arguments
    println!("Li0({}) = {}", z, z.li0());      // Li_0(z)
    println!("Li1({}) = {}", z, z.li1());      // Li_1(z)
    println!("Li2({}) = {}", z, z.li2());      // Li_2(z) (dilogarithm)
    println!("Li3({}) = {}", z, z.li3());      // Li_3(z) (trilogarithm)
    println!("Li4({}) = {}", z, z.li4());      // Li_4(z)
    println!("Li5({}) = {}", z, z.li5());      // Li_5(z)
    println!("Li6({}) = {}", z, z.li6());      // Li_6(z)
    println!("Li_{}({}) = {}", n, z, z.li(n)); // Li_n(z)
}
```


Notes
-----

The implementation of the real dilogarithm is an adaption of
[[arXiv:2201.01678](https://arxiv.org/abs/2201.01678)].

The implementation of the complex dilogarithm has been inspired by the
implementation in [SPheno](https://spheno.hepforge.org) and has been
translated to Rust.

The implementation of the general n-th order polylogarithm is an
adaption of [[arXiv:2010.09860](https://arxiv.org/abs/2010.09860)].


Citation
--------

~~~.bibtex
@software{polylog.rs,
    author  = {{Alexander Voigt}},
    title   = {{polylog.rs}},
    year    = {2022},
    version = {2.4.0},
    url     = {https://github.com/Expander/polylog},
    note    = {[License: LGPL-3.0-only]}
}
~~~


Copying
-------

Polylog is licenced under the GNU Lesser General Public License (GNU
LGPL) version 3.
