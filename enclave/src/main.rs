use std::net::{TcpStream};
use std::io::{Read, Result};
//use std::hash::{Hasher, BuildHasher};
use std::sync::{Arc, Mutex};
use std::ops::{DerefMut, Deref};
use byteorder::{ReadBytesExt, NetworkEndian};
use rand::thread_rng;
use rayon::scope;
use crossbeam::channel::{unbounded, bounded, Sender, Receiver};
use crate::cms::{CmsMap};
//use crate::cms::{CmsMap, CmsHasher, CmsBuildHasher};
use crate::parameters::*;
use crate::decryption::EncryptedReader;
use std::collections::HashMap;
use std::collections::hash_map::Entry;

mod chisq;
mod cms;
mod parameters;
mod decryption;

#[derive(Hash, Eq, PartialEq, Debug)]
struct Dat
{
    pub case_count: u16,
    pub control_count: u16,
}

impl Dat
{
	fn new(case_count: u16, control_count: u16) -> Dat
	{
		Dat
		{
			case_count: case_count,
			control_count: control_count
		}
    }
}

const ALLELE_HETEROZYGOUS: i8 = 1;
const ALLELE_HOMOZYGOUS: i8 = 2;
static mut STATIC_MAP: [[i16; WIDTH]; DEPTH] = [[0i16; WIDTH]; DEPTH];

/*
fn process_input_single_file_task(stream: &mut impl Read,  
                       txs: &[Sender<(Arc<Vec<u32>>, i8)>],
                       n_elems: usize, count: i8) {
    if cfg!(feature = "small-table"){
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
    } else {
        let mut buf = Vec::with_capacity(n_elems);
        for _ in 0..n_elems {
            let item = stream.read_u32::<NetworkEndian>().unwrap();
            buf.push(item);
        }
        let buf = Arc::new(buf);
        for tx in txs {
            tx.send((buf.clone(), count)).unwrap();
        }
    }
}
*/

fn process_input_single_file_task(stream: &mut impl Read, txs: &[Sender<Arc<Vec<u32>>>], n_elems: usize)
{
	let mut buf = Vec::with_capacity(n_elems);
	for _ in 0..n_elems
	{
		let item = stream.read_u32::<NetworkEndian>().unwrap();
		buf.push(item);
	}

	let buf = Arc::new(buf);
	for tx in txs
	{
		tx.send(buf.clone()).unwrap();
	}
}

/*
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
*/

fn process_input_files_task(stream: &mut impl Read,  
                            txs_meta: &[Sender<(usize, usize)>],
                            txs: &[Sender<Arc<Vec<u32>>>],
                            quota_txs: &[Receiver<()>],
                            n_files: usize,
             			    ) -> Result<()>
{
	for _ in 0..n_files
	{
		/*
		for quota_tx in quota_txs
		{
			quota_tx.recv().unwrap();
		}
		*/

		let patient_status = stream.read_u32::<NetworkEndian>()? as usize;
		let num_het_start = stream.read_u32::<NetworkEndian>()? as usize;
		let n_elems = stream.read_u32::<NetworkEndian>()? as usize;

		for tx in txs_meta
		{
			tx.send((patient_status, num_het_start)).unwrap();
		}

		process_input_single_file_task(stream, &txs[..], n_elems);
//		process_input_single_file_task(stream, &txs[..], n_elems);
	}

	/*
    for _ in 0..N_FILE_QUOTA
	{
		for quota_tx in quota_txs
		{
			quota_tx.recv().unwrap();
		}
	}
	*/
	Ok(())
}

fn process_tables_task(rx: &Receiver<(Arc<Vec<u32>>, i8)>,
                       map: &mut CmsMap,
                       n_elems: usize) {
    let rounds = {
        if cfg!(feature = "small-table"){
            (n_elems+BUF_SIZE-1)/BUF_SIZE
        } else {
            1
        }
    };

    for _ in 0..rounds {
        let (items, count) = rx.recv().unwrap();
        if cfg!(feature = "small-table")  {
            let positions = (*items).iter()
                .map(|x| map.cal_pos(*x) as u32)
                .collect::<Vec<_>>();
            for pos in positions {
                map.update_pos(pos, count as i16);
            }
        } else {
            let positions = (*items).iter()
                .map(|x| map.cal_pos(*x) as u32)
                .collect::<Vec<_>>();
            for i in 0..N_PARTITIONS {
                let range = (i as u32*PARTITION_SIZE as u32)..
                            ((i as u32+1)*PARTITION_SIZE as u32);
                for pos in &positions {
                    map.update_pos_in_range(*pos, count as i16, &range);
                }
            }
        }
    }

}

fn process_rhht_task(rx: &Receiver<Arc<Vec<u32>>>,
//                     map: &mut HashMap<u32, (u16, u16)>,
                     map: &mut HashMap<u32, Dat>,
                     patient_status: usize,
				     num_het_start: usize)
{
	let items = rx.recv().unwrap();
//	println!("1: {}\t{}", items.len(), num_het_start);

	for i in 0..num_het_start
	{
		if patient_status == 0
		{
			let mut ent = map.entry(items[i]).or_insert(Dat::new(0, 0));
			ent.control_count += 2;
//			*map.entry(items[i]).or_insert(Dat::new(0, 2)).control_count += 2;
//			(*map.get_mut(&items[i]).unwrap()).control_count += 2;
		}
		else
		{
			let mut ent = map.entry(items[i]).or_insert(Dat::new(0, 0));
			ent.case_count += 2;
			//*map.entry(items[i]).or_insert((2, 0)) += (2, 0);
		}
	}

//	println!("2: {}\t{}", items.len(), num_het_start);

	for i in num_het_start..items.len()
	{
		if patient_status == 0
		{
			let mut ent = map.entry(items[i]).or_insert(Dat::new(0, 0));
			ent.control_count += 1;
			//*map.entry(items[i]).or_insert((0, 1)) += (0, 1);
		}
		else
		{
			let mut ent = map.entry(items[i]).or_insert(Dat::new(0, 0));
			ent.case_count += 1;
			//*map.entry(items[i]).or_insert((1, 0)) += (1, 0);
		}
	}
}

fn main() {
    // connect to client
    let host = "localhost:1234";
    let mut stream = EncryptedReader::with_capacity(TCP_BUFFER_SIZE,
        TcpStream::connect(&host).expect("Tcp connect error."), &DUMMY_KEY);

    let n_files = stream.read_u32::<NetworkEndian>().unwrap() as usize;
    let stream = Arc::new(Mutex::new(stream));

    // initialize maps
//    let mut rng = thread_rng();
//    let maps = 
//        unsafe {
//            STATIC_MAP.iter_mut() 
//                .map(|m| Arc::new(Mutex::new(CmsMap::new(&mut rng, m))))
//                .collect::<Vec<_>>()
//        };

	// Test: Initialize HashMap
//	let mut rhht: HashMap<u32, (u16, u16)> = HashMap::new();
	let mut rhht = HashMap::new();
	let maps = vec![Arc::new(Mutex::new(rhht))];

    // initialize channels
    let (txs, rxs): (Vec<_>, Vec<_>) = (0..DEPTH)
        .map(|_| unbounded::<Arc<Vec<u32>>>())
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

/*
    for quota_tx in &quota_txs {
        for _ in 0..N_FILE_QUOTA {
            quota_tx.send(()).unwrap();
        }
    }
	*/

    let quota_txs = quota_txs.into_iter()
        .map(|x| Arc::new(Mutex::new(x)))
        .collect::<Vec<_>>();

    // set thread pool size
    rayon::ThreadPoolBuilder::new().num_threads(N_THREAD).build_global().unwrap();

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
/*
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
*/
		// Test: spawn map updaters
        for _ in 0..n_files {
            let iter = maps.iter()
                .zip(rxs.clone().into_iter())
                .zip(rxs_meta.clone().into_iter())
                .zip(quota_txs.clone().into_iter());

            for (((m, rx), rx_meta), quota_tx) in iter {
                s.spawn(move |_| {
                    let rx_meta = rx_meta.lock().unwrap();
                    let (patient_status, num_het_start) = rx_meta.recv().unwrap();
                    let rx = rx.lock().unwrap();
                    let mut m = m.lock().unwrap();

					// Test: Update HashMap
                    process_rhht_task(rx.deref(), m.deref_mut(), patient_status, num_het_start);

                    //let quota_tx = quota_tx.lock().unwrap();
                    //quota_tx.send(()).unwrap();
                });
            }
        }
    });

	for map in maps.into_iter()
	{
		let a =  Arc::try_unwrap(map).unwrap().into_inner().unwrap();
		for (key, value) in &a
		{
			let chisq = chisq::chi_sq(value.case_count, value.control_count, 2000, 2000);
			println!("{}\t{}\t{}\t{}", key, value.case_count, value.control_count, chisq);
		}
		//println!("{}", m);
		//println!("{}:\t{}\t{}", key, value.case_count, value.control_count);
	}

	// TODO: Second pass
	// TODO: RHHT/Priority Queue
}
