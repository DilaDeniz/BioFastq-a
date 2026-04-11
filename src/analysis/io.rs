use std::io::{self, BufWriter, Write};

use flate2::write::GzEncoder;
use flate2::Compression;

/// Open a file for writing; gzip-compress if the path ends in `.gz`.
pub fn open_writer(path: &str) -> io::Result<Box<dyn Write>> {
    let file = std::fs::File::create(path)?;
    if path.ends_with(".gz") {
        Ok(Box::new(BufWriter::new(GzEncoder::new(file, Compression::default()))))
    } else {
        Ok(Box::new(BufWriter::new(file)))
    }
}

/// Write one FASTQ record. For FASTA input (no quality), emits placeholder `I` scores.
pub fn write_fastq_record(
    writer: &mut dyn Write,
    id: &[u8],
    seq: &[u8],
    qual: Option<&[u8]>,
) -> io::Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    match qual {
        Some(q) => writer.write_all(q)?,
        None => writer.write_all(&vec![b'I'; seq.len()])?,
    }
    writer.write_all(b"\n")
}
