mod chisq;
mod csk;
mod cskf;
mod mcsk;
mod parameters;
mod decryption;
mod svd;
mod util;

use std::net::{TcpStream};
//use std::io::{Read, Result};
//use std::hash::{Hasher, BuildHasher};
use std::sync::{Arc, Mutex};
use std::ops::{DerefMut, Deref};
//use byteorder::{ReadBytesExt, NetworkEndian};
use rand::thread_rng;
use rayon::scope;
use crossbeam::channel::{unbounded, bounded, Sender, Receiver};
use crate::csk::{CskMap};
//use crate::cms::{CmsMap};
//use crate::cms::{CmsMap, CmsHasher, CmsBuildHasher};
use crate::parameters::*;
use crate::decryption::EncryptedReader;
use std::collections::HashMap;
use std::collections::BinaryHeap;

use std::env;
use std::fs::{read_dir, File};
use std::path::Path;
use std::net::{TcpListener};
use std::process::{Command, Child};
use std::io::{Write, Read, BufReader, Result, SeekFrom, Seek};
use std::mem::size_of;
use byteorder::{ReadBytesExt, WriteBytesExt, NetworkEndian, NativeEndian, LittleEndian, ByteOrder};

use std::cmp;
use std::cmp::Ordering;

use rand::seq::SliceRandom;

use ndarray::prelude::*;
use crate::svd::svdcomp;
use crate::util::matrix_ortho_proj1;


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

struct DatPcc
{
    pub ssqg: u16,
    pub dotprod: f32,
	pub sx: f32,
	pub pc_projections: Vec<f32>,
}

impl DatPcc
{
	fn new(ssqg: u16, dotprod: f32, sx: f32) -> DatPcc
	{
		DatPcc
		{
			ssqg: ssqg,
			dotprod: dotprod,
			sx: sx,
			pc_projections: Vec::new(),
		}
    }
}

#[derive(Hash, Eq, PartialEq, Ord, Debug)]
struct Snp(u32, u16);

impl PartialOrd for Snp
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering>
	{
        Some(other.1.cmp(&self.1))
    }
}

const ALLELE_HETEROZYGOUS: i16 = 1;
const ALLELE_HOMOZYGOUS: i16 = 2;
static mut STATIC_MAP: [[i16; WIDTH]; DEPTH] = [[0i16; WIDTH]; DEPTH];
//static mut STATIC_MAP: [[f32; WIDTH]; DEPTH] = [[0i16; WIDTH]; DEPTH];

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

fn process_sketch_task(rx: &Receiver<Arc<Vec<u32>>>,
                       map: &mut CskMap,
                       //map: &mut CmsMap,
                       patient_status: usize,
					   num_het_start: usize)
{
	let items = rx.recv().unwrap();

	let mut sign: i16 = 1;
	if patient_status == 0
	{
		sign = -1;
	}

	for i in 0..num_het_start
	{
		//map.update_pos(map.cal_pos(items[i]), ALLELE_HOMOZYGOUS * sign);
		map.update_pos_csk(items[i], ALLELE_HOMOZYGOUS * sign);
	}

	for i in num_het_start..items.len()
	{
		//map.update_pos(map.cal_pos(items[i]), ALLELE_HETEROZYGOUS * sign);
		map.update_pos_csk(items[i], ALLELE_HETEROZYGOUS * sign);
	}
}

fn process_rhht_task(rx: &Receiver<Arc<Vec<u32>>>,
                     map: &mut HashMap<u32, Dat>,
                     patient_status: usize,
				     num_het_start: usize)
{
	let items = rx.recv().unwrap();

	if patient_status == 0
	{
		for i in 0..num_het_start
		{
			let mut ent = map.entry(items[i]).or_insert(Dat::new(0, 0));
			ent.control_count += ALLELE_HOMOZYGOUS as u16;
		}

		for i in num_het_start..items.len()
		{
			let mut ent = map.entry(items[i]).or_insert(Dat::new(0, 0));
			ent.control_count += ALLELE_HETEROZYGOUS as u16;
		}
	}
	else
	{
		for i in 0..num_het_start
		{
			let mut ent = map.entry(items[i]).or_insert(Dat::new(0, 0));
			ent.case_count += ALLELE_HOMOZYGOUS as u16;
		}

		for i in num_het_start..items.len()
		{
			let mut ent = map.entry(items[i]).or_insert(Dat::new(0, 0));
			ent.case_count += ALLELE_HETEROZYGOUS as u16;
		}
	}
}

fn process_rhht_task_temp(items: Vec<u32>, map: &mut HashMap<u32, Dat>)
{
	let patient_status = items[0];
	let num_het_start = items[1];

	if patient_status == 0
	{
		for i in 2..(num_het_start + 2) as usize
		{
			if map.contains_key(&items[i as usize])
			{
				if let Some(x) = map.get_mut(&items[i as usize])
				{
					x.control_count += ALLELE_HOMOZYGOUS as u16;
				}
			}
		}

		for i in (num_het_start + 2)..items.len() as u32
		{
			if map.contains_key(&items[i as usize])
			{
				if let Some(x) = map.get_mut(&items[i as usize])
				{
					x.control_count += ALLELE_HETEROZYGOUS as u16;
				}
			}
		}
	}
	else
	{
		for i  in 2..(num_het_start + 2) as usize
		{
			if map.contains_key(&items[i as usize])
			{
				if let Some(x) = map.get_mut(&items[i as usize])
				{
					x.case_count += ALLELE_HOMOZYGOUS as u16;
				}
			}
		}

		for i in (num_het_start + 2)..items.len() as u32
		{
			if map.contains_key(&items[i as usize])
			{
				if let Some(x) = map.get_mut(&items[i as usize])
				{
					x.case_count += ALLELE_HETEROZYGOUS as u16;
				}
			}
		}
	}
}

fn process_mcsk_task_temp(items: Vec<u32>, mcsk: &mut mcsk::Mcsk, phenotypes: &mut Array<f32, Ix1>, file_idx: usize)
{
	let patient_status = items[0];
	let num_het_start = items[1];

	phenotypes[file_idx] = -1.0;
	if patient_status == 0
	{
		phenotypes[file_idx] = 1.0;
	}

	for i in 2..(num_het_start + 2)
	{
		mcsk.mcsk_update(items[i as usize], file_idx, 2.0);
	}

	for i in (num_het_start + 2)..items.len() as u32
	{
		mcsk.mcsk_update(items[i as usize], file_idx, 1.0);
	}
}

fn main()
{
    // Connect to client
    let host = "localhost:1234";
    let mut stream = EncryptedReader::with_capacity(TCP_BUFFER_SIZE,
        TcpStream::connect(&host).expect("Tcp connect error."), &DUMMY_KEY);

    let n_files = stream.read_u32::<NetworkEndian>().unwrap() as usize;
    let stream = Arc::new(Mutex::new(stream));

	// Initialize sketches
    let mut rng = thread_rng();
    let maps = 
        unsafe {
            STATIC_MAP.iter_mut() 
                //.map(|m| Arc::new(Mutex::new(CmsMap::new(&mut rng, m))))
                .map(|m| Arc::new(Mutex::new(CskMap::new(&mut rng, m))))
                .collect::<Vec<_>>()
        };

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
        // Spawn input processor
        s.spawn(move |_| {
            let mut stream = stream.lock().unwrap();
            process_input_files_task(stream.deref_mut(), 
                                     &txs_meta[..], 
                                     &txs[..], 
                                     &quota_rxs[..],
                                     n_files).unwrap();
        });

        // Spawn sketch updaters
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

                    process_sketch_task(rx.deref(), m.deref_mut(), patient_status, num_het_start);

                    //let quota_tx = quota_tx.lock().unwrap();
                    //quota_tx.send(()).unwrap();
                });
            }
        }

		// Update RHHT
		/*
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

					// Update RHHT
                    process_rhht_task(rx.deref(), m.deref_mut(), patient_status, num_het_start);

                    //let quota_tx = quota_tx.lock().unwrap();
                    //quota_tx.send(()).unwrap();
                });
            }
        }
		*/
    });

/*
	let mut f = File::open("/home/ckockan/chr1_uniq/chr1_uniq.ckz0").unwrap();
	let mut buffer: Vec<u8> = Vec::new();
	f.read_to_end(&mut buffer).unwrap();
	
	let mut snps: Vec<u32> = Vec::new();
	for i in 0..(buffer.len() / 4)
	{
		snps.push(0);
	}
	LittleEndian::read_u32_into(&buffer, &mut snps);

	// Initialize Priority Queue as Heap
	let mut heap: BinaryHeap<Snp> = BinaryHeap::new();

	let maps = maps.into_iter().map(|m| Arc::try_unwrap(m).ok().unwrap().into_inner().unwrap()).collect::<Vec<_>>();

	for snp in &snps
	{
		let mut values = Vec::new();
		for m in maps.iter()
		{
			values.push(m.csk_get_median_odd(*snp as u64));
		}

		// Sort values
		values.sort();

		// Get median of values
		let mut median: i16 = 0;
		if DEPTH % 2 == 1
		{
			median = values[DEPTH / 2];
		}
		else
		{
			if values[DEPTH / 2] < -200
			{
				median = values[DEPTH / 2 - 1];
			}
			else if values[DEPTH / 2 - 1] > 200
			{
				median = values[DEPTH / 2];
			}
			else
			{
				median = (values[DEPTH / 2 - 1] + values[DEPTH / 2]) / 2;
			}
		}

		if heap.len() >= 10000
		{
			let mut val = heap.peek_mut().unwrap();
			if val.1 < median.abs() as u16
			{
				*val = Snp(*snp, median.abs() as u16);
			}
		}
		else
		{
			heap.push(Snp(*snp, median.abs() as u16));
		}
	}

	/* DEBUG
	println!("Heap size: {}", heap.len());
	for _ in 0..9990
	{
		heap.pop();
	}

	for _ in 0..10
	{
		println!("{:?}", heap.pop());
	}
	*/

	// Initialize RHHT
//	let mut rhht = HashMap::new();
//	for x in heap.iter()
//	{
//		rhht.insert(x.0, Dat::new(0, 0));
//	}
*/
    let mut rng = thread_rng();
	let mut mcsk = mcsk::Mcsk::new(&mut rng);
	let mut phenotypes: Array<f32, Ix1> = Array::zeros(2000);
	let mut file_idx: usize = 0;

	let dir_path = "/home/ckockan/chr1/";
    let dir = Path::new(dir_path);
    if dir.is_dir()
	{
		println!("Second pass ...");

		let mut rng = rand::thread_rng();
		let mut file_names = read_dir(dir).unwrap().collect::<Vec<_>>();
		file_names.as_mut_slice().shuffle(&mut rng);

		for (_i, entry) in file_names.into_iter().enumerate()
		{  
            println!("File {}...", _i);
            let path = entry.unwrap().path();
			let mut f = File::open(path).unwrap();
			let mut buffer: Vec<u8> = Vec::new();
			f.read_to_end(&mut buffer).unwrap();
	
			let mut input: Vec<u32> = Vec::new();
			for i in 0..(buffer.len() / 4)
			{
				input.push(0);
			}
			LittleEndian::read_u32_into(&buffer, &mut input);

			//process_rhht_task_temp(input, &mut rhht);
			process_mcsk_task_temp(input, &mut mcsk, &mut phenotypes, file_idx);
			file_idx = file_idx + 1;
        }
    }
	else
	{
		panic!("Empty directory.");
    }

	println!("{:?}", mcsk.mcsk);
	//mcsk.mcsk_print();
	mcsk.mcsk_mean_centering();
	println!();
	//mcsk.mcsk_print();
	println!("{:?}", mcsk.mcsk);

	// DEBUG: svd
	let mut u: Array<f32, Ix1> = Array::ones(MCSK_WIDTH);

	let mut A = mcsk.mcsk.view_mut();

	if MCSK_WIDTH < MCSK_DEPTH
	{
		let mut S: Array<f32, Ix1> = Array::zeros(MCSK_WIDTH);

		let mut Q: Array<f32, Ix2> = Array::zeros((MCSK_WIDTH, MCSK_WIDTH));

		// Compute SVD A = USV^T; V stored in Q
		//let mut retval: i32 = svdcomp_t(A, MCSK_DEPTH, MCSK_WIDTH, S, Q);
		svdcomp(A.slice_mut(s![.., ..A.ncols() - 1]), S.view_mut(), Q.view_mut());
		Q = Q.reversed_axes();

		// Copy MCSK_NUM_PC rows of Q to A.
		A.slice_mut(s![0..MCSK_NUM_PC, ..A.ncols() - 1]).assign(&Q.slice(s![0..MCSK_NUM_PC, ..]));
//		for i in 0..MCSK_NUM_PC
//		{
//			for j in 0..MCSK_WIDTH
//			{
//				A[[i, j]] = Q[[i, j]];
//			}
//		}
		
		// Compute VV^T * phenotype vector and VV^T * all one vector
		for i in 0..MCSK_WIDTH
		{
			A[[MCSK_NUM_PC, i]] = 0.0;
		}

		matrix_ortho_proj1(A.view_mut(), phenotypes.view(), MCSK_NUM_PC, MCSK_WIDTH);

		for i in 0..MCSK_WIDTH
		{
			// Should be replaced by daxpy
			A[[MCSK_NUM_PC, i]] = phenotypes[i] - A[[MCSK_NUM_PC, i]];
		}

		for i in 0..MCSK_WIDTH
		{
			phenotypes[i] = A[[MCSK_NUM_PC, i]];
		}

		for i in 0..MCSK_WIDTH
		{
			A[[MCSK_NUM_PC + 1, i]] = 0.0;
		}

		matrix_ortho_proj1(A.view_mut(), u.view(),  MCSK_NUM_PC, MCSK_WIDTH);

		for i in 0..MCSK_WIDTH
		{
			u[i] = A[[MCSK_NUM_PC + 1, i]];
		}
	}

	// Keep only the first k rows of V, now stroed in the first k rows of A
	let mut enclave_eig: Array<f32, Ix2> = Array::zeros((MCSK_NUM_PC, MCSK_WIDTH));
	for pc in 0..MCSK_NUM_PC
	{
		for i in 0..MCSK_WIDTH
		{
			enclave_eig[[pc, i]] = A[[pc, i]];
		}
	}
	println!("{:?}", enclave_eig);
	// END DEBUG: svd	
/*
	let mut diff_vec = Vec::new();
	for (key, value) in &rhht
	{
		let chisq = chisq::chi_sq(value.case_count, value.control_count, 2000, 2000);
		//println!("{}\t{}\t{}\t{}", key, value.case_count, value.control_count, chisq);
		let abs_diff = (value.case_count as i32 - value.control_count as i32).abs();
		diff_vec.push(abs_diff);
	}

	diff_vec.sort();
	diff_vec.reverse();
	for x in diff_vec.iter()
	{
		println!("{}", x);
	}
	println!("{}", diff_vec.len());
*/
	/*
	for map in maps.into_iter()
	{
		let a =  Arc::try_unwrap(map).unwrap().into_inner().unwrap();
		for (key, value) in &a
		{
			let chisq = chisq::chi_sq(value.case_count, value.control_count, 2000, 2000);
			println!("{}\t{}\t{}\t{}", key, value.case_count, value.control_count, chisq);
		}
	}
	*/
}
