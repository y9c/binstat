use bam::record::Aux;
use rust_htslib::faidx;
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::PathBuf;
use std::str::FromStr;

fn count_base(c: u8) -> char {
    // A = 65, a = 97
    // C = 67, c = 99
    // G = 71, g = 103
    // T = 84, t = 116
    // U = 85, u = 117
    match c {
        65 => 'A',
        97 => 'A',
        85 => 'T',
        117 => 'T',
        67 => 'C',
        99 => 'C',
        71 => 'G',
        103 => 'G',
        84 => 'T',
        116 => 'T',
        _ => 'N',
    }
}

pub fn run(bam_path: PathBuf, fasta_path: PathBuf, region_path: PathBuf, win: usize, step: usize) {
    //let args: Vec<String> = env::args().collect();
    //let bam_path = &args[1];
    //let fasta_path = &args[2];
    //let region_path: String;
    //if args.len() >= 4 {
    //    region_path = args[3].clone();
    //} else {
    //    region_path = fasta_path.clone() + ".fai";
    //}
    //let win = 100;
    //let step = 10000;

    let file = File::open(&region_path).unwrap();
    let reader = BufReader::new(file);

    let dna_bases = &vec!['A', 'C', 'G', 'T'];
    let mut gc_background: HashMap<usize, usize> = HashMap::new();
    let mut yf_counter: HashMap<usize, usize> = HashMap::new();
    let mut zf_counter: HashMap<usize, usize> = HashMap::new();
    let fa_reader = &faidx::Reader::from_path(&fasta_path).unwrap();
    let mut bam = bam::IndexedReader::from_path(&bam_path).unwrap();
    bam.set_reference(&fasta_path).expect("Set ref");

    for line in reader.lines() {
        match line {
            Ok(line) => {
                let mut region = line.split('\t');
                let chrom: &str;
                let start: usize;
                let end: usize;
                if region_path.ends_with("bed") {
                    chrom = region.next().unwrap();
                    start = usize::from_str(region.next().unwrap()).unwrap();
                    end = usize::from_str(region.next().unwrap()).unwrap();
                } else {
                    chrom = region.next().unwrap();
                    start = 0;
                    end = usize::from_str(region.next().unwrap()).unwrap();
                }
                // let chrom = "X";
                // let start = 202908;
                // let end = 5502908;
                //let end = i32::from_str(region.next().unwrap()).unwrap();
                eprintln!("{}:{}-{}", chrom, start, end);

                let mut span_start = start;
                let mut span_end = start + win;
                let mut span_index = 0;
                while span_start < end {
                    span_index += 1;
                    if span_index % 10000 == 0 {
                        eprintln!(
                            "Finished {}M regions on genome",
                            span_index as f64 / (1000000.0 / 100.0 / 1000.0)
                        );
                    }
                    let fa_string = fa_reader
                        .fetch_seq(&chrom, span_start as usize, (span_end - 1) as usize)
                        .unwrap();
                    let base_counter = dna_bases
                        .clone()
                        .into_iter()
                        .map(|b| fa_string.iter().filter(|&c| count_base(*c) == b).count())
                        .collect::<Vec<_>>();
                    let total = base_counter.iter().sum::<usize>();
                    let gc = (base_counter[1] + base_counter[2]) as f64 / total as f64;
                    let gci = (gc * 1000.0).round() as usize;
                    gc_background
                        .entry(gci)
                        .and_modify(|e| *e += 1)
                        .or_insert(1);
                    if total > 0 {
                        bam.fetch((chrom, span_start as i64, span_end as i64))
                            .unwrap();
                        for record in bam.records() {
                            let r = record.unwrap();
                            let yf = match r.aux(b"Yf") {
                                Ok(value) => {
                                    let v = if let Aux::I8(v) = value { v } else { todo!() };
                                    Some(v)
                                }
                                Err(_e) => None,
                            };
                            let zf = match r.aux(b"Zf") {
                                Ok(value) => {
                                    let v = if let Aux::I8(v) = value { v } else { todo!() };
                                    Some(v)
                                }
                                Err(_e) => None,
                            };
                            if yf.is_some() & zf.is_some() {
                                let yf_count = yf.unwrap() as usize;
                                let zf_count = zf.unwrap() as usize;
                                yf_counter
                                    .entry(gci)
                                    .and_modify(|e| *e += yf_count)
                                    .or_insert(yf_count);
                                zf_counter
                                    .entry(gci)
                                    .and_modify(|e| *e += zf_count)
                                    .or_insert(zf_count);
                            }
                        }
                    }
                    span_start += step;
                    span_end = span_start + win;
                }
            }
            Err(..) => {}
        }
    }
    for (k, v) in yf_counter {
        println!(
            "{}\t{}\t{}\t{}",
            k as f64 / 1000.0,
            v,
            zf_counter[&k],
            gc_background[&k]
        );
    }
}
