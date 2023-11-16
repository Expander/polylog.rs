use num::complex::Complex;
use std::fs::File;
use std::io::{BufReader, BufRead, Error, ErrorKind};
use std::path::PathBuf;


#[macro_export]
macro_rules! assert_eq_float {
    ($a:expr, $b:expr, $eps:expr) => {
        if ($a - $b).abs() >= $eps*(1.0 + $a.abs().max($b.abs())) {
            println!("Numbers differ by more than {}: {} != {}", $eps, $a, $b);
            assert!(false);
        }
    }
}


#[macro_export]
macro_rules! assert_eq_complex {
    ($a:expr, $b:expr, $eps:expr) => {
        assert_eq_float!($a.re, $b.re, $eps);
        assert_eq_float!($a.im, $b.im, $eps);
    }
}


pub trait CLn<T> {
    fn cln(&self) -> T;
}


impl CLn<Complex<f64>> for Complex<f64> {
    fn cln(&self) -> Complex<f64> {
        Complex::new(
            if self.re == 0. { 0. } else { self.re },
            if self.im == 0. { 0. } else { self.im },
        ).ln()
    }
}


fn data_path(filename: &str) -> PathBuf {
    let mut path = PathBuf::from(file!());
    path.pop();
    path.pop();
    path.pop();
    path.push("tests");
    path.push("data");
    path.push(filename);
    path
}


pub fn read_data_file(filename: &str) -> Result<Vec<(Complex<f64>, Complex<f64>)>, std::io::Error> {
    let file = File::open(data_path(filename))?;
    let br = BufReader::new(file);
    let mut vec = Vec::new();

    for line in br.lines() {
        let vals = line?
            .split_whitespace()
            .map(|s| s.parse::<f64>().unwrap())
            .collect::<Vec<f64>>();

        if vals.len() != 4 {
            return Err(Error::new(ErrorKind::UnexpectedEof, "line does not contain 4 real numbers"));
        }

        vec.push((Complex::new(vals[0], vals[1]), Complex::new(vals[2], vals[3])));
    }

    Ok(vec)
}
