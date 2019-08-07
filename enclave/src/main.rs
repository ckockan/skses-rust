use std::net::{TcpStream};
use std::io::{Read, Result, BufReader};
use std::hash::{Hasher, BuildHasher};
use byteorder::{ReadBytesExt, NetworkEndian, NativeEndian};
use rand::thread_rng;
use rayon::prelude::*;
use crate::cms::{CmsMap, CmsHasher, CmsBuildHasher};

mod cms;

const ALLELE_HETEROZYGOUS: i16 = 1;
const ALLELE_HOMOZYGOUS: i16 = 2;

const WIDTH: usize = 0x100000/4; // 0.25 MB
//const WIDTH: usize = 0x200000; // 2 MB
const DEPTH: usize = 8;
const N_THREAD: usize = 4;

const PARTITION: bool = false;
const N_PARTITIONS: usize = 8;
const PARTITION_LEN: usize = WIDTH/N_PARTITIONS;
fn process_single_file(stream: &mut impl Read, 
                       maps: &mut [CmsMap]) -> Result<()> {
    let patient_status = stream.read_u32::<NativeEndian>()?;
    let num_het_start = stream.read_u32::<NetworkEndian>()? as usize;
    let n_elems = stream.read_u32::<NetworkEndian>()?;
    let mut file_buf = Vec::with_capacity(n_elems as usize);
    for _ in 0..n_elems {
        let rs_id_uint = stream.read_u32::<NetworkEndian>()?;
        file_buf.push(rs_id_uint);
    }
    let count_hom = match patient_status==0 {
        true => -ALLELE_HOMOZYGOUS,
        false => ALLELE_HOMOZYGOUS,
    }; 
    let count_het = match patient_status==0 {
        true => -ALLELE_HETEROZYGOUS,
        false => ALLELE_HETEROZYGOUS,
    }; 
    let _ = maps.par_iter_mut()
        .map(|m| {
            for item in &file_buf[..num_het_start] {
                m.update(*item, count_het);
            }
            for item in &file_buf[num_het_start..] {
                m.update(*item, count_hom);
            }
        })
    .collect::<Vec<_>>();
    Ok(())
}

fn main() {
    rayon::ThreadPoolBuilder::new().num_threads(N_THREAD).build_global().unwrap();
    let host = "localhost:1234";
    let mut stream = BufReader::new(
        TcpStream::connect(&host).expect("Tcp connect error."));
    let n_files = stream.read_u32::<NetworkEndian>().unwrap();
    let mut rng = thread_rng();
    let mut maps = (0..DEPTH)
        .map(|_| CmsMap::new(WIDTH, &mut rng))
        .collect::<Vec<_>>();
    for _ in 0..n_files {
        process_single_file(&mut stream, &mut maps[..])
            .expect("Error recieving files.");
    }

}
