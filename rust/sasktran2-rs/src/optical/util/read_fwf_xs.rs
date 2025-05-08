use ndarray::*;
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader};
use std::ops::AddAssign;
use std::path::PathBuf;

/// Raw storage structure that holds the cross section data from something like the
/// HITRAN fwf files, or the LBLRTM files.  This is basically lists of cross section data as a function
/// of wavenumber at a given temperature and pressure.  We assume wavenumber is specified in
/// vacuum cm^-1, and the cross section is in cm^2/molecule.
pub struct XsecListAtConditions {
    pub xsec: Vec<Array1<f64>>,
    pub wvnum: Vec<Array1<f64>>,
    pub params: HashMap<String, Vec<f64>>,
}

/// Let's us concatenate two XsecDatabase objects together.  This is useful for combining
/// databases from different pressure/temperatures together
impl AddAssign for XsecListAtConditions {
    fn add_assign(&mut self, other: Self) {
        self.xsec.extend(other.xsec);
        self.wvnum.extend(other.wvnum);

        for (param, vals) in other.params {
            self.params.entry(param).or_default().extend(vals);
        }
    }
}

#[allow(dead_code)]
struct Header {
    pub short_mol_name: String,
    pub wvnum_start: f64,
    pub wvnum_end: f64,
    pub num_points: i64,
    pub temperature: f64,
    pub zero: f64,
    pub wvnum_space: f64,
    pub pressure: f64,
}

impl Header {
    fn new() -> Self {
        Header {
            short_mol_name: String::new(),
            wvnum_start: 0.0,
            wvnum_end: 0.0,
            num_points: 0,
            temperature: 0.0,
            zero: 0.0,
            wvnum_space: 0.0,
            pressure: 0.0,
        }
    }
}

fn read_hitran_header(line: &str) -> Header {
    let short_mol_name = line[10..20].trim().to_string();
    let wvnum_start: f64 = line[20..31].trim().parse().unwrap();
    let wvnum_end: f64 = line[31..42].trim().parse().unwrap();
    let num_points: i64 = line[42..50].trim().parse().unwrap();
    let temperature: f64 = line[50..58].trim().parse().unwrap();
    let zero: f64 = line[59..62].trim().parse().unwrap();
    // let wvnum_space: f64 = line[62..71].trim().parse().unwrap();  // unsure on this one...
    let pressure: f64 = line[71..80].trim().parse().unwrap();

    let wvnum_space = (wvnum_end - wvnum_start) / ((num_points - 1) as f64);

    Header {
        short_mol_name,
        wvnum_start,
        wvnum_end,
        num_points,
        temperature,
        zero,
        wvnum_space,
        pressure,
    }
}

pub fn read_fwf_xsec(path: PathBuf) -> Option<XsecListAtConditions> {
    let file = fs::File::open(path).unwrap();
    let reader = BufReader::new(file);

    let mut all_xs: Vec<Array1<f64>> = vec![];
    let mut all_wvnum: Vec<Array1<f64>> = vec![];

    let mut xs: Vec<f64> = vec![];
    let mut wvnum: Vec<f64> = vec![];
    let mut pressure: Vec<f64> = vec![];
    let mut temperature: Vec<f64> = vec![];

    let mut header = Header::new();
    let mut cur_wvnum = 0.0;
    let mut cur_index = 0;

    for line_result in reader.lines() {
        let line = line_result.ok()?;

        // Check if this is header information
        if line.len() == 102 {
            header = read_hitran_header(&line);
            cur_wvnum = header.wvnum_start;
        } else {
            for i in 0..10 {
                let start = i * 10;
                let end = start + 10;
                let value: f64 = line[start..end].trim().parse().unwrap();
                xs.push(value);
                wvnum.push(cur_wvnum);
                cur_wvnum += header.wvnum_space;

                cur_index += 1;
                if cur_index >= header.num_points {
                    all_xs.push(Array1::from(xs.clone()));
                    all_wvnum.push(Array1::from(wvnum.clone()));
                    pressure.push(header.pressure);
                    temperature.push(header.temperature);

                    xs.clear();
                    wvnum.clear();
                    break;
                }
            }
        }
    }

    let mut params = HashMap::new();
    params.insert("pressure".to_string(), pressure);
    params.insert("temperature".to_string(), temperature);

    Some(XsecListAtConditions {
        xsec: all_xs,
        wvnum: all_wvnum,
        params,
    })
}

pub fn read_fwf_folder(folder: PathBuf) -> Option<XsecListAtConditions> {
    let paths = fs::read_dir(folder).unwrap();

    let mut combined_ds = XsecListAtConditions {
        xsec: vec![],
        wvnum: vec![],
        params: HashMap::new(),
    };

    for path in paths {
        let path = path.unwrap().path();
        if path.extension().unwrap() == "xsc" {
            let dbase = read_fwf_xsec(path.clone());

            combined_ds += dbase?;
        }
    }

    Some(combined_ds)
}
