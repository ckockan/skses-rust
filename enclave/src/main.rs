use std::net::{TcpStream};
use std::io::{Read, Result, BufReader};
use std::hash::{Hasher, BuildHasher};
use std::sync::{Arc, Mutex};
use std::ops::{DerefMut, Deref};
use byteorder::{ReadBytesExt, NetworkEndian};
use rand::thread_rng;
use rayon::prelude::*;
use crossbeam::channel::{unbounded, bounded, Sender, Receiver};
use crate::cms::{CmsMap, CmsHasher, CmsBuildHasher};

mod cms;

const ALLELE_HETEROZYGOUS: i8 = 1;
const ALLELE_HOMOZYGOUS: i8 = 2;

const WIDTH: usize = 0x100000/4; // 0.25 MB
//const WIDTH: usize = 0x200000; // 2 MB
const DEPTH: usize = 8;
const N_THREAD: usize = 5;
const N_FILE_QUOTA: usize = 6;
const BUF_SIZE: usize = 2000; 

fn process_input_single_file_task(stream: &mut impl Read,  
                       txs: &[Sender<(Arc<Vec<u32>>, i8)>],
                       n_elems: usize, count: i8) {

    let rounds = (n_elems+BUF_SIZE-1)/BUF_SIZE;
    for i in 0..rounds {
        let mut buf = Vec::with_capacity(BUF_SIZE);
        for _ in (i*BUF_SIZE)..usize::min(n_elems, (i+1)*BUF_SIZE){
            let item = stream.read_u32::<NetworkEndian>().unwrap();
            buf.push(item);
        }
        let buf = Arc::new(buf);
        for tx in txs {
            tx.send((buf.clone(), count)).unwrap();
        }
    }
}

fn process_input_files_task(stream: &mut impl Read,  
                      txs_meta: &[Sender<(usize, usize)>],
                      txs: &[Sender<(Arc<Vec<u32>>, i8)>],
                      quota_tx: &Receiver<()>,
                      n_files: usize,
                      ) -> Result<()> {

    for _ in 0..n_files {
        //println!("len: {}", quota_tx.len());
        quota_tx.recv().unwrap();
        let patient_status = stream.read_u32::<NetworkEndian>()?;
        let num_het_start = stream.read_u32::<NetworkEndian>()? as usize;
        let n_elems = stream.read_u32::<NetworkEndian>()? as usize;
        let n_hom = num_het_start;
        let n_het = n_elems - num_het_start;

        for tx in txs_meta {
            tx.send((n_hom, n_het)).unwrap();
        }

        let count_hom = match patient_status==0 {
            true => -ALLELE_HOMOZYGOUS,
            false => ALLELE_HOMOZYGOUS,
        }; 
        let count_het = match patient_status==0 {
            true => -ALLELE_HETEROZYGOUS,
            false => ALLELE_HETEROZYGOUS,
        }; 

        process_input_single_file_task(stream, &txs[..], n_hom, count_hom);
        process_input_single_file_task(stream, &txs[..], n_het, count_het);
    }
    for _ in 0..N_FILE_QUOTA {
        quota_tx.recv().unwrap();
    }
    Ok(())
}

fn process_tables_task(rx: &Receiver<(Arc<Vec<u32>>, i8)>,
                       map: &mut CmsMap,
                       n_elems: usize) {
    let rounds = (n_elems+BUF_SIZE-1)/BUF_SIZE;
    for _ in 0..rounds {
        let (items, count) = rx.recv().unwrap();
        for item in &*items {
            map.update(*item, count as i16);
        }
    }

}

fn main() {
    // set thread pool size
    rayon::ThreadPoolBuilder::new().num_threads(N_THREAD).build_global().unwrap();

    // connect to client
    let host = "localhost:1234";
    let mut stream = BufReader::with_capacity(0x10000,
        TcpStream::connect(&host).expect("Tcp connect error."));
    let n_files = stream.read_u32::<NetworkEndian>().unwrap() as usize;
    let stream = Arc::new(Mutex::new(stream));

    // initialize maps
    let mut rng = thread_rng();
    let mut maps = (0..DEPTH)
        .map(|_| CmsMap::new(WIDTH, &mut rng))
        .collect::<Vec<_>>();

    // initialize channels
    let (txs, rxs): (Vec<_>, Vec<_>) = (0..DEPTH)
        .map(|_| unbounded::<(Arc<Vec<u32>>, i8)>())
        .unzip();

    let (txs_meta, rxs_meta): (Vec<_>, Vec<_>) = (0..DEPTH)
        .map(|_| unbounded::<(usize, usize)>())
        .unzip();

    // create quota to limit the number of files to read in at a time
    let (quota_tx, quota_rx) = bounded::<()>(N_FILE_QUOTA);
    for _ in 0..N_FILE_QUOTA {
        quota_tx.send(()).unwrap();
    }

    // spawn input processers
    rayon::scope(|s| {
        s.spawn(move |_| {
            let quota_rx = quota_rx;
            let mut stream = stream.lock().unwrap();
            process_input_files_task(stream.deref_mut(), 
                                     &txs_meta[..], 
                                     &txs[..], 
                                     &quota_rx,
                                     n_files).unwrap();
        });

        let rxs = rxs.into_iter()
            .map(|x| Arc::new(Mutex::new(x)))
            .collect::<Vec<_>>();

        let rxs_meta = rxs_meta.into_iter()
            .map(|x| Arc::new(Mutex::new(x)))
            .collect::<Vec<_>>();
        for _ in 0..n_files {
            let _ = maps.par_iter_mut()
                .zip_eq(rxs.clone().into_par_iter())
                .zip_eq(rxs_meta.clone().into_par_iter())
                .map(|((m, rx), rx_meta)| {
                    let rx_meta = rx_meta.lock().unwrap();
                    let (n_hom, n_het) = rx_meta.recv().unwrap();

                    let rx = rx.lock().unwrap();
                    process_tables_task(rx.deref(), m, n_hom);
                    process_tables_task(rx.deref(), m, n_het);
                })
            .collect::<Vec<_>>();
            quota_tx.send(()).unwrap();
        }


    });

}
