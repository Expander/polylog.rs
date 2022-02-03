//! Polylog
//! =======
//!
//! The Polylog package provides Rust implementations of real and
//! complex polylogarithms, including the dilogarithm and
//! trilogarithm.
//!
//! # Example:
//! ```
//! extern crate num;
//! extern crate polylog;
//! use num::complex::Complex;
//! use polylog::{Li2, Li3, Li4, Li5, Li6};
//!
//! fn main() {
//!     let x = 1.0;
//!     let z = Complex::new(1.0, 1.0);
//!     println!("Li2({}) = {}", x, x.li2()); // Re[Li_2(x)] (real dilogarithm)
//!     println!("Li2({}) = {}", z, z.li2()); // Li_2(z)     (complex dilogarithm)
//!     println!("Li3({}) = {}", x, x.li3()); // Re[Li_3(x)] (real trilogarithm)
//!     println!("Li3({}) = {}", z, z.li3()); // Li_3(z)     (complex trilogarithm)
//!     println!("Li4({}) = {}", x, x.li4()); // Re[Li_4(x)]
//!     println!("Li4({}) = {}", z, z.li4()); // Li_4(z)
//!     println!("Li5({}) = {}", z, z.li5()); // Li_5(z)
//!     println!("Li6({}) = {}", z, z.li6()); // Li_6(z)
//! }
//! ```


extern crate num;

mod cln;
mod digamma;
mod li0;
mod li1;
mod li2;
mod li3;
mod li4;
mod li5;
mod li6;
mod li;
mod zeta;

pub use self::li0::Li0;
pub use self::li1::Li1;
pub use self::li2::Li2;
pub use self::li3::Li3;
pub use self::li4::Li4;
pub use self::li5::Li5;
pub use self::li6::Li6;
pub use self::li::Li;
