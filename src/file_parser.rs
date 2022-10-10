use bio::io::fasta;
use bio::io::fastq;
use std::fs;
use std::fs::metadata;

/// get_genomes(genomes_path: &String)
/// --
/// obtains a list of genome objects given a genome directory,
/// or a genome file (i.e. fasta or fastq)
///
/// Parameters
/// ----------
///
/// Returns
/// -------
///
/// Notes
/// -----
///
/// Examples
/// --------
///
pub fn get_genomes(genomes_path: &String) {
    let md = metadata(genomes_path).unwrap();
    if md.is_dir() {
        let paths = fs::read_dir(genomes_path).unwrap();
        for path in paths {
            //println!("{:?}", path.unwrap().path().to_str().unwrap());
            parse_genome_file(&path.unwrap().path().to_str().unwrap().to_string());
        }
    } else {
        parse_genome_file(genomes_path);
    }
}

/// parse_genome_file(genomes_file_path: &String)
/// --
/// obtains a list of genome objects given a genome file (i.e. fasta or fastq)
///
/// Parameters
/// ----------
///
/// Returns
/// -------
///
/// Notes
/// -----
///
/// Examples
/// --------
///
fn parse_genome_file(genomes_file_path: &String) {
    println!("\nParsing file: {}\n", genomes_file_path);
    if genomes_file_path.ends_with(".fa")
        || genomes_file_path.ends_with(".fasta")
        || genomes_file_path.ends_with(".fna")
    {
        fasta_parser(&genomes_file_path);
    } else if genomes_file_path.ends_with(".fq") || genomes_file_path.ends_with(".fastq") {
        fastq_parser(&genomes_file_path);
    }
}

/// fasta_parser(file_path: &String)
/// --
/// parses a fasta file.
///
/// Parameters
/// ----------
///
/// Returns
/// -------
///
/// Notes
/// -----
///
/// Examples
/// --------
///
fn fasta_parser(file_path: &String) {
    let mut nb_reads = 0;
    let mut nb_bases = 0;
    let reader = fasta::Reader::from_file(file_path).unwrap();
    //let fa_reader = reader.from_file("fake_fasta.fa");
    for result in reader.records() {
        let result_data = &result.unwrap();

        nb_reads += 1;
        nb_bases += result_data.seq().len();
        println!("{:?}", result_data.id());
        //println!("Genome Length: {:?}", result_data.seq());
    }
    println!("Number of reads: {}", nb_reads);
    println!("Number of bases: {}", nb_bases);
}

/// fastq_parser(file_path: &String)
/// --
/// parses a fastq file.
///
/// Parameters
/// ----------
///
/// Returns
/// -------
///
/// Notes
/// -----
///
/// Examples
/// --------
///
fn fastq_parser(file_path: &String) {
    let mut nb_reads = 0;
    let mut nb_bases = 0;
    let reader = fastq::Reader::from_file(file_path).unwrap();
    for result in reader.records() {
        let record = &result.unwrap();
        nb_reads += 1;
        nb_bases += record.seq().len();
        //println!("{:?}", result_data.id());
        //println!("Genome Length: {:?}", record.seq());
    }
    println!("Number of reads: {}", nb_reads);
    println!("Number of bases: {}", nb_bases);
}
