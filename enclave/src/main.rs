mod chisq;
mod csk;
mod cskf;
mod decryption;
mod mcsk;
mod parameters;
mod svd;
mod util;

use std::env;
use std::fs::{read_dir, File};
use std::path::Path;
use std::net::{TcpListener, TcpStream};
use std::process::{Command, Child};
use std::io::{Write, Read, BufRead, BufReader, Result, SeekFrom, Seek};
use std::mem::size_of;
use std::sync::{Arc, Mutex};
use std::ops::{DerefMut, Deref};
use std::collections::HashMap;
use std::collections::BinaryHeap;
use std::cmp;
use std::cmp::Ordering;
use rayon::scope;
use rand::thread_rng;
use rand::seq::SliceRandom;
use byteorder::{ReadBytesExt, WriteBytesExt, NetworkEndian, NativeEndian, LittleEndian, ByteOrder};
use crossbeam::channel::{unbounded, bounded, Sender, Receiver};
use ndarray::prelude::*;
use ndarray::{Array1, Array2, ArrayView1, ArrayView2, ArrayViewMut1, ArrayViewMut2};
use crate::parameters::*;
use crate::csk::{CskMap};
use crate::cskf::{CskfMap};
use crate::decryption::EncryptedReader;
use crate::svd::svdcomp;
use crate::util::*;

pub struct Parameters
{
	pub host: String,
	pub csk_width: usize,
	pub csk_depth: usize,
	pub mcsk_width: usize,
	pub mcsk_depth: usize,
	pub mcsk_num_pc: usize,
}

impl Parameters
{
	fn new(host: String, csk_width: usize, csk_depth: usize, mcsk_width: usize, mcsk_depth: usize, mcsk_num_pc: usize) -> Parameters
	{
		Parameters
		{
			host: host,
			csk_width: csk_width,
			csk_depth: csk_depth,
			mcsk_width: mcsk_width,
			mcsk_depth: mcsk_depth,
			mcsk_num_pc: mcsk_num_pc,
		}
	}
}

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
	fn new(ssqg: u16, dotprod: f32, sx: f32, num_pc: usize) -> DatPcc
	{
		DatPcc
		{
			ssqg: ssqg,
			dotprod: dotprod,
			sx: sx,
			pc_projections: vec![0.0; num_pc],
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

fn process_query_file_task(stream: &mut impl Read, snps: &mut Vec<u32>) -> Result<()>
{
	let num_snps_uniq = stream.read_u32::<NetworkEndian>()? as usize;
	for _ in 0..num_snps_uniq
	{
		let item = stream.read_u32::<NetworkEndian>()?;
		snps.push(item);
	}
	Ok(())
}

fn process_sketch_task(rx: &Receiver<Arc<Vec<u32>>>,
                       map: &mut CskMap,
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
		map.update_pos_csk(items[i], ALLELE_HOMOZYGOUS * sign);
	}

	for i in num_het_start..items.len()
	{
		map.update_pos_csk(items[i], ALLELE_HETEROZYGOUS * sign);
	}
}

fn process_cskf_task(rx: &Receiver<Arc<Vec<u32>>>,
                       map: &mut CskfMap,
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
		map.update_pos_cskf(items[i], (ALLELE_HOMOZYGOUS as f32) * (sign as f32));
	}

	for i in num_het_start..items.len()
	{
		map.update_pos_cskf(items[i], (ALLELE_HETEROZYGOUS as f32) * (sign as f32));
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

fn process_rhht_pcc_task(items: Vec<u32>, map: &mut HashMap<u32, DatPcc>, num_pc: usize, phenotypes: ArrayView1<f32>, enclave_eig: ArrayView2<f32>, u: ArrayView1<f32>, file_idx: usize)
{
	let patient_status = items[0];
	let num_het_start = items[1];

	for i in 2..(num_het_start + 2)
	{
		let mut ent = map.entry(items[i as usize]).or_insert(DatPcc::new(0, 0.0, 0.0, num_pc));
		ent.ssqg += 4;
		ent.dotprod += phenotypes[file_idx] * 2.0;
		ent.sx += (2.0 * (1.0 - u[file_idx]));
		for pc in 0..ent.pc_projections.len()
		{
			ent.pc_projections[pc] += enclave_eig[[pc, file_idx]] * 2.0;
		}
	}

	for i in (num_het_start + 2)..items.len() as u32
	{
		let mut ent = map.entry(items[i as usize]).or_insert(DatPcc::new(0, 0.0, 0.0, num_pc));
		ent.ssqg += 1;
		ent.dotprod += phenotypes[file_idx];
		ent.sx += 1.0 - u[file_idx];
		for pc in 0..ent.pc_projections.len()
		{
			ent.pc_projections[pc] += enclave_eig[[pc, file_idx]];
		}
	}
}

fn process_rhht_pcc_task_temp(items: Vec<u32>, map: &mut HashMap<u32, DatPcc>, num_pc: usize, phenotypes: ArrayView1<f32>, enclave_eig: ArrayView2<f32>, u: ArrayView1<f32>, file_idx: usize)
{
	let patient_status = items[0];
	let num_het_start = items[1];

	for i in 2..(num_het_start + 2) as usize
	{
		if map.contains_key(&items[i as usize])
		{
			if let Some(x) = map.get_mut(&items[i as usize])
			{
				x.ssqg += 4;
				x.dotprod += phenotypes[file_idx] * 2.0;
				x.sx += 2.0 * (1.0 - u[file_idx]);
				for pc in 0..x.pc_projections.len()
				{
					x.pc_projections[pc] += enclave_eig[[pc, file_idx]] * 2.0;
				}
			}
		}
	}

	for i in (num_het_start + 2)..items.len() as u32
	{
		if map.contains_key(&items[i as usize])
		{
			if let Some(x) = map.get_mut(&items[i as usize])
			{
				x.ssqg += 1;
				x.dotprod += phenotypes[file_idx];
				x.sx += 1.0 - u[file_idx];
				for pc in 0..x.pc_projections.len()
				{
					x.pc_projections[pc] += enclave_eig[[pc, file_idx]];
				}
			}
		}
	}
}

fn process_mcsk_task_temp(items: Vec<u32>, mcsk: &mut mcsk::Mcsk, phenotypes: &mut Array1<f32>, file_idx: usize)
{
	let patient_status = items[0];
	let num_het_start = items[1];

	phenotypes[file_idx] = 0.0;
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

fn process_mcsk_task(rx: &Receiver<Arc<Vec<u32>>>, mcsk: &mut mcsk::Mcsk, patient_status: usize, num_het_start: usize, phenotypes: &mut Array1<f32>, file_idx: usize)
{
	let items = rx.recv().unwrap();
	
	phenotypes[file_idx] = 0.0;
	if patient_status == 0
	{
		phenotypes[file_idx] = 1.0;
	}

	for i in 0..num_het_start
	{
		mcsk.mcsk_update(items[i as usize], file_idx, 2.0);
	}

	for i in num_het_start..items.len()
	{
		mcsk.mcsk_update(items[i as usize], file_idx, 1.0);
	}
}

fn process_rhht_pcc_task_strm(rx: &Receiver<Arc<Vec<u32>>>, map: &mut HashMap<u32, DatPcc>, num_het_start: usize, num_pc: usize, phenotypes: &mut Array1<f32>, enclave_eig: &mut Array2<f32>, u: &mut Array1<f32>, file_idx: usize)
{
	let items = rx.recv().unwrap();

	for i in 0..num_het_start as usize
	{
		if map.contains_key(&items[i as usize])
		{
			if let Some(x) = map.get_mut(&items[i as usize])
			{
				x.ssqg += 4;
				x.dotprod += phenotypes[file_idx] * 2.0;
				x.sx += 2.0 * (1.0 - u[file_idx]);
				for pc in 0..x.pc_projections.len()
				{
					x.pc_projections[pc] += enclave_eig[[pc, file_idx]] * 2.0;
				}
			}
		}
	}

	for i in num_het_start..items.len()
	{
		if map.contains_key(&items[i as usize])
		{
			if let Some(x) = map.get_mut(&items[i as usize])
			{
				x.ssqg += 1;
				x.dotprod += phenotypes[file_idx];
				x.sx += 1.0 - u[file_idx];
				for pc in 0..x.pc_projections.len()
				{
					x.pc_projections[pc] += enclave_eig[[pc, file_idx]];
				}
			}
		}
	}
}

fn process_query_task(rx: &Receiver<Arc<Vec<u32>>>) -> Vec<u32>
{
	let mut v: Vec<u32> = Vec::new();

	let items = rx.recv().unwrap();

	for i in 0..items.len()
	{
		v.push(items[i]);
	}
	return v;
}

fn main()
{
	// Parse configuration file
	// TODO: Don't hardcode the config file path, make it a command line parameter
	let mut params = Parameters::new(String::new(), 0, 0, 0, 0, 0);
	let file_config = File::open("/home/ckockan/skses-rust/config.txt").unwrap();
	let reader_config = BufReader::new(file_config);
	for line in reader_config.lines()
	{
		let ln = line.unwrap();
		let tokens: Vec<&str> = ln.split('=').collect();

		match tokens[0]
		{
			"HOST" => params.host = tokens[1].to_string(),
			"CSK_WIDTH" => params.csk_width = tokens[1].parse().unwrap(),
			"CSK_DEPTH" => params.csk_depth = tokens[1].parse().unwrap(),
			"MCSK_WIDTH" => params.mcsk_width = tokens[1].parse().unwrap(),
			"MCSK_DEPTH" => params.mcsk_depth = tokens[1].parse().unwrap(),
			"MCSK_NUM_PC" => params.mcsk_num_pc = tokens[1].parse().unwrap(),
			_ => println!("Warning: Unrecognized Parameter"),
		}
	}

    // Connect to client
	let mut stream = EncryptedReader::with_capacity(TCP_BUFFER_SIZE, TcpStream::connect(&params.host).expect("Tcp connect error."), &DUMMY_KEY);

	// Receive the number of files to be received from the service provider
    let n_files = stream.read_u32::<NetworkEndian>().unwrap() as usize;

	// Setup the stream
    let stream = Arc::new(Mutex::new(stream));

	// Initialize sketches
    let mut rng = thread_rng();
    //let csk_maps = (0..params.csk_depth).map(|_| Arc::new(Mutex::new(CskMap::new(&mut rng, &params)))).collect::<Vec<_>>();
    let cskf_maps = (0..params.csk_depth).map(|_| Arc::new(Mutex::new(CskfMap::new(&mut rng, &params)))).collect::<Vec<_>>();

    // Initialize channels
    let (txs, rxs): (Vec<_>, Vec<_>) = (0..params.csk_depth)
        .map(|_| unbounded::<Arc<Vec<u32>>>())
        .unzip();

    let rxs = rxs.into_iter()
        .map(|x| Arc::new(Mutex::new(x)))
        .collect::<Vec<_>>();

    let (txs_meta, rxs_meta): (Vec<_>, Vec<_>) = (0..params.csk_depth)
        .map(|_| unbounded::<(usize, usize)>())
        .unzip();

    let rxs_meta = rxs_meta.into_iter()
        .map(|x| Arc::new(Mutex::new(x)))
        .collect::<Vec<_>>();

    // create quota to limit the number of files to read in at a time
    let (quota_txs, quota_rxs): (Vec<_>, Vec<_>) = (0..params.csk_depth)
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

	let mcsk = Arc::new(Mutex::new(mcsk::Mcsk::new(&mut rng, &params)));
	let mut phenotypes: Array1<f32> = Array1::zeros(2000);
	let mut file_idx: usize = 0;
	let snps: Arc<Mutex<Vec<u32>>> = Arc::new(Mutex::new(Vec::new()));
	let mut enclave_eig: Array2<f32> = Array2::zeros((params.mcsk_num_pc, params.mcsk_width));
	let mut u: Array1<f32> = Array1::ones(params.mcsk_width);
	let mut rhht_pcc = HashMap::new();

    scope(|s|
	{
        // Spawn input processor
		let snps_input_processor = snps.clone();
		let (t_sync, r_sync) = unbounded::<()>();
        s.spawn(move |_|
		{
			let mut stream = stream.lock().unwrap();
			process_input_files_task(stream.deref_mut(), &txs_meta[0..1], &txs[0..1], &quota_rxs[0..1], n_files).unwrap();

    		stream.deref_mut().read_u32::<NetworkEndian>().unwrap() as usize;
            process_input_files_task(stream.deref_mut(), &txs_meta[..], &txs[..], &quota_rxs[..], n_files).unwrap();

			{
				let mut snps = snps_input_processor.lock().unwrap();
				process_query_file_task(stream.deref_mut(), snps.deref_mut()).unwrap();
			}
			t_sync.send(()).unwrap();

    		stream.deref_mut().read_u32::<NetworkEndian>().unwrap() as usize;
            process_input_files_task(stream.deref_mut(), &txs_meta[0..1], &txs[0..1], &quota_rxs[0..1], n_files).unwrap();
        });

		println!("First Pass: Updating MCSK ...");
		for _ in 0..n_files
		{
			let rx_meta = rxs_meta[0].lock().unwrap();
			let (patient_status, num_het_start) = rx_meta.recv().unwrap();
			let rx = rxs[0].lock().unwrap();
			let mut m = mcsk.lock().unwrap();
			process_mcsk_task(rx.deref(), m.deref_mut(), patient_status, num_het_start, &mut phenotypes, file_idx);
			file_idx = file_idx + 1;
		}

		//let mcsk = Arc::try_unwrap(mcsk).unwrap().into_inner().unwrap();
		let mut mcsk = mcsk.lock().unwrap();

		// Perform mean centering prior to SVD
		mcsk.mcsk_mean_centering();

		// SVD
		println!("Performing SVD ...");
		let mut a = mcsk.mcsk.view_mut();
		if params.mcsk_width < params.mcsk_depth
		{
			let mut s: Array1<f32> = Array1::zeros(params.mcsk_width);
			let mut q: Array2<f32> = Array2::zeros((params.mcsk_width, params.mcsk_width));

			// Compute SVD A = U * S * V^T
			// V is stored in Q
			svdcomp(a.slice_mut(s![.., ..a.ncols() - 1]), s.view_mut(), q.view_mut());

			// Copy MCSK_NUM_PC rows of Q to A.
			a.slice_mut(s![0..params.mcsk_num_pc, ..a.ncols() - 1]).assign(&q.slice(s![0..params.mcsk_num_pc, ..]));

			// Compute V * V^T * phenotype_vector and V * V^T * all_ones_vector
			for i in 0..params.mcsk_width
			{
				a[[params.mcsk_num_pc, i]] = 0.0;
			}

			matrix_ortho_proj1(a.view_mut(), phenotypes.view(), params.mcsk_num_pc, params.mcsk_width);

			for i in 0..params.mcsk_width
			{
				// Should be replaced by daxpy (?, ask Kaiyuan)
				a[[params.mcsk_num_pc, i]] = phenotypes[i] - a[[params.mcsk_num_pc, i]];
			}

			for i in 0..params.mcsk_width
			{
				phenotypes[i] = a[[params.mcsk_num_pc, i]];
			}

			for i in 0..params.mcsk_width
			{
				a[[params.mcsk_num_pc, i]] = 0.0;
			}

			matrix_ortho_proj1(a.view_mut(), u.view(), params.mcsk_num_pc, params.mcsk_width);

			for i in 0..params.mcsk_width
			{
				u[i] = a[[params.mcsk_num_pc, i]];
			}
		}

		// Keep only the first k rows of V, now stroed in the first k rows of A
		for pc in 0..params.mcsk_num_pc
		{
			for i in 0..params.mcsk_width
			{
				enclave_eig[[pc, i]] = a[[pc, i]];
			}
		}

        // Spawn sketch updaters
		// scope to make sure everything within scope|r| finishes before proceeding
		println!("Second Pass: Updating CSKF ...");
    	scope(|r|
		{
			for _ in 0..n_files
			{
				let iter = cskf_maps.iter().zip(rxs.clone().into_iter()).zip(rxs_meta.clone().into_iter()).zip(quota_txs.clone().into_iter());
				for (((m, rx), rx_meta), quota_tx) in iter
				{
					r.spawn(move |_|
					{
						let rx_meta = rx_meta.lock().unwrap();
						let (patient_status, num_het_start) = rx_meta.recv().unwrap();
						let rx = rx.lock().unwrap();
						let mut m = m.lock().unwrap();

						//process_sketch_task(rx.deref(), m.deref_mut(), patient_status, num_het_start);
						process_cskf_task(rx.deref(), m.deref_mut(), patient_status, num_het_start);

						//let quota_tx = quota_tx.lock().unwrap();
						//quota_tx.send(()).unwrap();
					});
				}
			}

		});

		// Initialize Priority Queue as Heap
		println!("Populating Heap ...");
		let mut heap: BinaryHeap<Snp> = BinaryHeap::new();

		let maps = cskf_maps.into_iter().map(|m| Arc::try_unwrap(m).ok().unwrap().into_inner().unwrap()).collect::<Vec<_>>();

		r_sync.recv().unwrap();

		let mut snps = snps.lock().unwrap();
		for i in 0..snps.len()
		{
			let mut values = Vec::new();
			for m in maps.iter()
			{
				values.push(m.cskf_get_median_odd(snps[i] as u64));
			}

			// Sort values
			values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

			let mut median: f32 = 0.0;
			if params.csk_depth % 2 == 1
			{
				median = values[params.csk_depth / 2];
			}
			else
			{
				if values[params.csk_depth / 2] < -200.0
				{
					median = values[params.csk_depth / 2 - 1];
				}
				else if values[params.csk_depth / 2 - 1] > 200.0
				{
					median = values[params.csk_depth / 2];
				}
				else
				{
					median = (values[params.csk_depth / 2 - 1] + values[params.csk_depth / 2]) / 2.0;
				}
			}

			if heap.len() >= 10000
			{
				let mut val = heap.peek_mut().unwrap();
				if val.1 < median.abs() as u16
				{
					*val = Snp(snps[i], median.abs() as u16);
				}
			}
			else
			{
				heap.push(Snp(snps[i], median.abs() as u16));
			}
		}

		println!("Initializing RHHT-PCC ...");
		for x in heap.iter()
		{
			rhht_pcc.insert(x.0, DatPcc::new(0, 0.0, 0.0, params.mcsk_num_pc));
		}

		println!("Updating RHHT-PCC ...");
		file_idx = 0;
		for _ in 0..n_files
		{
			let rx_meta = rxs_meta[0].lock().unwrap();
			let (patient_status, num_het_start) = rx_meta.recv().unwrap();
			let rx = rxs[0].lock().unwrap();
			process_rhht_pcc_task_strm(rx.deref(), &mut rhht_pcc, num_het_start, params.mcsk_num_pc, &mut phenotypes, &mut enclave_eig, &mut u, file_idx);
			file_idx = file_idx + 1;
		}

		println!("Reporting results ...");
		let mut sy: f32 = 0.0;
		let mut sy2: f32 = 0.0;
		for x in phenotypes.iter()
		{
			sy = sy + x;
			sy2 = sy2 + (x * x);
		}

		for (key, value) in &rhht_pcc
		{
			let cat_chisq = chisq::cat_chi_sq(value.ssqg, 2000, value.dotprod, value.sx, sy, sy2, &value.pc_projections);
			println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", key, cat_chisq, value.ssqg, value.dotprod, value.sx, sy, sy2, &value.pc_projections[0], &value.pc_projections[1]);
		}
	});
}
