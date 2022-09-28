use rust_htslib::faidx;
use rust_htslib::{bam, bam::Read};

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

fn main() {
    let fasta_path = "./test/test.fa";
    let bam_path = "./test/test.bam";
    let win = 1000000;
    let dna_bases = &vec!['A', 'C', 'G', 'T'];

    let fa_reader = &faidx::Reader::from_path(&fasta_path).unwrap();
    let mut bam = bam::IndexedReader::from_path(&bam_path).unwrap();

    let chrom = "X";
    let start = 202908;
    let end = start + win;

    let fa_string = fa_reader
        .fetch_seq(&chrom, start as usize, (end - 1) as usize)
        .unwrap();
    let base_counter = dna_bases
        .clone()
        .into_iter()
        .map(|b| fa_string.iter().filter(|&c| count_base(*c) == b).count())
        .collect::<Vec<_>>();
    let gc = (base_counter[1] + base_counter[2]) as f64 / base_counter.iter().sum::<usize>() as f64;
    println!("{}:{}-{} \n {:?}", chrom, start, start + win, gc);

    bam.fetch((chrom, start as i64, end as i64)).unwrap();
    //println!("{}", bam.records().count());
    let mut c_fwd = 0;
    let mut c_rev = 0;
    for read in bam.records() {
        let read = read.unwrap();
        match read.aux(b"YZ") {
            Ok(n) => {
                // 43 = '+', 45 = '-'
                match n {
                    bam::record::Aux::Char(43) => c_fwd += 1,
                    bam::record::Aux::Char(45) => c_rev += 1,
                    _ => (),
                }
            }
            Err(..) => {
                //println!("no NH tag")
            }
        }
    }
    println!(
        "{}:{}-{} \n {:?} {:?}",
        chrom,
        start,
        start + win,
        c_fwd,
        c_rev
    );
}
