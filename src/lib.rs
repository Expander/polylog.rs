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
//! use polylog::li2::Li2;
//! use polylog::li3::Li3;
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

pub mod li2;
pub mod li3;
