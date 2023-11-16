//! Polylog
//! =======
//!
//! The Polylog package provides Rust implementations of real and
//! complex polylogarithms, including the dilogarithm and
//! trilogarithm.
//!
//! # Example:
//! ```
//! use num::complex::Complex;
//! use polylog::{Li, Li0, Li1, Li2, Li3, Li4, Li5, Li6};
//!
//! fn main() {
//!     let x = 1.0;
//!     let z = Complex::new(1.0, 1.0);
//!     let n = 10;
//!
//!     // real polylogarithms for real arguments
//!     println!("Li0({}) = {}", x, x.li0());      // Re[Li_0(x)]
//!     println!("Li1({}) = {}", x, x.li1());      // Re[Li_1(x)]
//!     println!("Li2({}) = {}", x, x.li2());      // Re[Li_2(x)] (dilogarithm)
//!     println!("Li3({}) = {}", x, x.li3());      // Re[Li_3(x)] (trilogarithm)
//!     println!("Li4({}) = {}", x, x.li4());      // Re[Li_4(x)]
//!     println!("Li_{}({}) = {}", n, x, x.li(n)); // Re[Li_n(x)]
//!
//!     // complex polylogarithms for complex arguments
//!     println!("Li0({}) = {}", z, z.li0());      // Li_0(z)
//!     println!("Li1({}) = {}", z, z.li1());      // Li_1(z)
//!     println!("Li2({}) = {}", z, z.li2());      // Li_2(z) (dilogarithm)
//!     println!("Li3({}) = {}", z, z.li3());      // Li_3(z) (trilogarithm)
//!     println!("Li4({}) = {}", z, z.li4());      // Li_4(z)
//!     println!("Li5({}) = {}", z, z.li5());      // Li_5(z)
//!     println!("Li6({}) = {}", z, z.li6());      // Li_6(z)
//!     println!("Li_{}({}) = {}", n, z, z.li(n)); // Li_n(z)
//! }
//! ```


mod cln;
mod li0;
mod li1;
mod li2;
mod li3;
mod li4;
mod li5;
mod li6;
mod li;

pub use self::li0::Li0;
pub use self::li1::Li1;
pub use self::li2::Li2;
pub use self::li3::Li3;
pub use self::li4::Li4;
pub use self::li5::Li5;
pub use self::li6::Li6;
pub use self::li::Li;
