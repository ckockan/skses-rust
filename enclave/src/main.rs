use std::net::{TcpStream};
use std::io::{Read, Result, BufReader};
//use std::hash::{Hasher, BuildHasher};
use std::sync::{Arc, Mutex};
use std::ops::{DerefMut, Deref};
use byteorder::{ReadBytesExt, NetworkEndian};
use rand::thread_rng;
use rayon::scope;
use crossbeam::channel::{unbounded, bounded, Sender, Receiver};
use crate::cms::{CmsMap};
//use crate::cms::{CmsMap, CmsHasher, CmsBuildHasher};

mod cms;

const ALLELE_HETEROZYGOUS: i8 = 1;
const ALLELE_HOMOZYGOUS: i8 = 2;

const WHOLE_FILE: bool = false;
const WIDTH: usize = 0x40000;  // 0.25 MB

//const WHOLE_FILE: bool = true;
//const WIDTH: usize = 0x200000; // 2 MB

const N_THREAD: usize = 3;

const BUF_SIZE: usize = 2000; 
const DEPTH: usize = 8;
const N_FILE_QUOTA: usize = 6;
const PARTITION_SIZE: usize = 0x40000; 
const N_PARTITIONS: usize = (WIDTH+PARTITION_SIZE-1)/PARTITION_SIZE;

fn process_input_single_file_task(stream: &mut impl Read,  
                       txs: &[Sender<(Arc<Vec<u32>>, i8)>],
                       n_elems: usize, count: i8) {
    if WHOLE_FILE {
        let mut buf = Vec::with_capacity(n_elems);
        for _ in 0..n_elems {
            let item = stream.read_u32::<NetworkEndian>().unwrap();
            buf.push(item);
        }
        let buf = Arc::new(buf);
        for tx in txs {
            tx.send((buf.clone(), count)).unwrap();
        }

    } else {
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
}

fn process_input_files_task(stream: &mut impl Read,  
                      txs_meta: &[Sender<(usize, usize)>],
                      txs: &[Sender<(Arc<Vec<u32>>, i8)>],
                      quota_txs: &[Receiver<()>],
                      n_files: usize,
                      ) -> Result<()> {
    for _ in 0..n_files {
        for quota_tx in quota_txs {
            quota_tx.recv().unwrap();
        }
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
        for quota_tx in quota_txs {
            quota_tx.recv().unwrap();
        }
    }
    Ok(())
}

fn process_tables_task(rx: &Receiver<(Arc<Vec<u32>>, i8)>,
                       map: &mut CmsMap,
                       n_elems: usize) {
    let rounds = {
        if WHOLE_FILE {
            1
        } else {
            (n_elems+BUF_SIZE-1)/BUF_SIZE
        }
    };

    for _ in 0..rounds {
        let (items, count) = rx.recv().unwrap();
        if WIDTH <= 0x40000  {
            let positions = (*items).iter()
                .map(|x| map.cal_pos(*x) as u32)
                .collect::<Vec<_>>();
            for pos in positions {
                map.update_pos(pos as usize, count as i16);
            }
        } else {
            let positions = (*items).iter()
                .map(|x| map.cal_pos(*x) as u32)
                .collect::<Vec<_>>();
            for i in 0..N_PARTITIONS {
                let range = (i*PARTITION_SIZE)..((i+1)*PARTITION_SIZE);
                for pos in &positions {
                    map.update_pos_in_range(*pos as usize, count as i16, &range);
                }
            }
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
    let maps = (0..DEPTH)
        .map(|_| Arc::new(Mutex::new(CmsMap::new(WIDTH, &mut rng))))
        .collect::<Vec<_>>();

    // initialize channels
    let (txs, rxs): (Vec<_>, Vec<_>) = (0..DEPTH)
        .map(|_| unbounded::<(Arc<Vec<u32>>, i8)>())
        .unzip();

    let rxs = rxs.into_iter()
        .map(|x| Arc::new(Mutex::new(x)))
        .collect::<Vec<_>>();

    let (txs_meta, rxs_meta): (Vec<_>, Vec<_>) = (0..DEPTH)
        .map(|_| unbounded::<(usize, usize)>())
        .unzip();

    let rxs_meta = rxs_meta.into_iter()
        .map(|x| Arc::new(Mutex::new(x)))
        .collect::<Vec<_>>();

    // create quota to limit the number of files to read in at a time
    let (quota_txs, quota_rxs): (Vec<_>, Vec<_>) = (0..DEPTH)
        .map(|_| bounded::<()>(N_FILE_QUOTA))
        .unzip();

    for quota_tx in &quota_txs {
        for _ in 0..N_FILE_QUOTA {
            quota_tx.send(()).unwrap();
        }
    }
    let quota_txs = quota_txs.into_iter()
        .map(|x| Arc::new(Mutex::new(x)))
        .collect::<Vec<_>>();


    scope(|s| {
        // spawn input processor
        s.spawn(move |_| {
            let mut stream = stream.lock().unwrap();
            process_input_files_task(stream.deref_mut(), 
                                     &txs_meta[..], 
                                     &txs[..], 
                                     &quota_rxs[..],
                                     n_files).unwrap();
        });

        // spawn table updaters
        for _ in 0..n_files {
            let iter = maps.iter()
                .zip(rxs.clone().into_iter())
                .zip(rxs_meta.clone().into_iter())
                .zip(quota_txs.clone().into_iter());

            for (((m, rx), rx_meta), quota_tx) in iter {
                s.spawn(move |_| {
                    let rx_meta = rx_meta.lock().unwrap();
                    let (n_hom, n_het) = rx_meta.recv().unwrap();
                    let rx = rx.lock().unwrap();
                    let mut m = m.lock().unwrap();
                    process_tables_task(rx.deref(), m.deref_mut(), n_hom);
                    process_tables_task(rx.deref(), m.deref_mut(), n_het);
                    let quota_tx = quota_tx.lock().unwrap();
                    quota_tx.send(()).unwrap();
                });
            }
        }
    });
}
