use std::net::{TcpStream};
use std::io::{Read, Result, BufReader};
use std::hash::{Hasher, BuildHasher};
use std::sync::{mpsc::{channel, Receiver, Sender}, Arc, Mutex};
use std::ops::DerefMut;
use byteorder::{ReadBytesExt, NetworkEndian};
use rand::thread_rng;
use rayon::prelude::*;
use crate::cms::{CmsMap, CmsHasher, CmsBuildHasher};

mod cms;

const ALLELE_HETEROZYGOUS: i8 = 1;
const ALLELE_HOMOZYGOUS: i8 = 2;

const WIDTH: usize = 0x100000/4; // 0.25 MB
//const WIDTH: usize = 0x200000; // 2 MB
const DEPTH: usize = 8;
const N_THREAD: usize = 5;
const BUF_SIZE: usize = 2000; 

fn process_inputs_task(stream: &mut impl Read,  
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

fn process_single_file<T>(stream: Arc<Mutex<T>>, 
                       maps: &mut [CmsMap]) -> Result<()>
where T: Read + Send {
    let (patient_status, num_het_start, n_elems) = {
        let mut stream = stream.lock().unwrap();
        let patient_status = stream.read_u32::<NetworkEndian>()?;
        let num_het_start = stream.read_u32::<NetworkEndian>()? as usize;
        let n_elems = stream.read_u32::<NetworkEndian>()? as usize;
        (patient_status, num_het_start, n_elems)
    };
    let n_hom = num_het_start;
    let n_het = n_elems - num_het_start;
    let count_hom = match patient_status==0 {
        true => -ALLELE_HOMOZYGOUS,
        false => ALLELE_HOMOZYGOUS,
    }; 
    let count_het = match patient_status==0 {
        true => -ALLELE_HETEROZYGOUS,
        false => ALLELE_HETEROZYGOUS,
    }; 
    let (txs, rxs): (Vec<_>, Vec<_>) = (0..DEPTH)
        .map(|_| channel::<(Arc<Vec<u32>>, i8)>())
        .unzip();
    
    rayon::scope(|s| {
        s.spawn(move |_| {
            let mut stream = stream.lock().unwrap();
            let stream_ref = stream.deref_mut();
            process_inputs_task(stream_ref, &txs[..], n_hom, count_hom);
            process_inputs_task(stream_ref, &txs[..], n_het, count_het);
        });
    });

    let _ = maps.par_iter_mut()
        .zip_eq(rxs.into_par_iter())
        .map(|(m, rx)| {
            process_tables_task(&rx, m, n_hom);
            process_tables_task(&rx, m, n_het);
        })
    .collect::<Vec<_>>();
    Ok(())
}

fn main() {
    rayon::ThreadPoolBuilder::new().num_threads(N_THREAD).build_global().unwrap();

    let host = "localhost:1234";
    let mut stream = BufReader::with_capacity(0x10000,
        TcpStream::connect(&host).expect("Tcp connect error."));
    let n_files = stream.read_u32::<NetworkEndian>().unwrap();
    let mut rng = thread_rng();
    let mut maps = (0..DEPTH)
        .map(|_| CmsMap::new(WIDTH, &mut rng))
        .collect::<Vec<_>>();
    let stream = Arc::new(Mutex::new(stream));
    for _ in 0..n_files {
        process_single_file(stream.clone(), &mut maps[..])
            .expect("Error recieving files.");
    }

}
