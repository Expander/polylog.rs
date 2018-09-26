//! Polylog
//! =======
//!
//! The Polylog package provides Rust implementations of real and
//! complex polylogarithms.
//!
//! # Example:
//! ```
//! extern crate num;
//! extern crate polylog;
//! use num::complex::Complex;
//! use polylog::{Li2, Li3};
//!
//! fn main() {
//!     let x = 1.0;
//!     let z = Complex::new(1.0, 1.0);
//!     println!("Li2({}) = {}", x, x.li2());
//!     println!("Li2({}) = {}", z, z.li2());
//!     println!("Li3({}) = {}", z, z.li3());
//! }
//! ```


extern crate num;

mod li2;
mod li3;

pub use self::li2::Li2;
pub use self::li3::Li3;
