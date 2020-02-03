extern crate num;
use num::complex::Complex;
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::path::PathBuf;

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
        let line = line?;
        let vals = line
            .trim()
            .split("\t")
            .collect::<Vec<&str>>()
            .iter()
            .map(|s| s.parse::<f64>().unwrap())
            .collect::<Vec<f64>>();

        assert!(vals.len() == 4);

        vec.push((Complex::new(vals[0], vals[1]), Complex::new(vals[2], vals[3])));
    }

    Ok(vec)
}
