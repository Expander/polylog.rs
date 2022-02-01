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
use polylog::{Li2, Li3, Li4, Li5, Li6};

fn main() {
    let x = 1.0;
    let z = Complex::new(1.0, 1.0);
    println!("Li2({}) = {}", x, x.li2());
    println!("Li2({}) = {}", z, z.li2());
    println!("Li3({}) = {}", x, x.li3());
    println!("Li3({}) = {}", z, z.li3());
    println!("Li4({}) = {}", z, z.li4());
    println!("Li5({}) = {}", z, z.li5());
    println!("Li6({}) = {}", z, z.li6());
}
```


Notes
-----

The implementation of the complex dilogarithm has been inspired by the
implementation in [SPheno](https://spheno.hepforge.org) and has been
translated to Rust.


Copying
-------

Polylog is licenced under the GNU Lesser General Public License (GNU
LGPL) version 3.
