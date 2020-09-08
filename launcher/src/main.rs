use std::env;
use std::fs::{read_dir, File};
use std::path::Path;
use std::net::{TcpListener};
use std::process::{Command, Child};
use std::io::{Write, Read, BufReader, Result, SeekFrom, Seek};
use std::mem::size_of;
use byteorder::{ReadBytesExt, WriteBytesExt, NetworkEndian, NativeEndian};
use crate::parameters::*;
use crate::encryption::EncryptedWriter;
use rand::seq::SliceRandom;

mod encryption;
mod parameters;

fn send_single_file(file: &mut impl Read, stream: &mut impl Write, n_elems: usize) -> Result<()>
{
	let patient_status = file.read_u32::<NativeEndian>()?;
	stream.write_u32::<NetworkEndian>(patient_status)?;

	let num_het_start = file.read_u32::<NativeEndian>()?;
	stream.write_u32::<NetworkEndian>(num_het_start)?;

	stream.write_u32::<NetworkEndian>(n_elems as u32)?;

	for _ in 0..n_elems
	{
		let rs_id_uint = file.read_u32::<NativeEndian>()?;
		stream.write_u32::<NetworkEndian>(rs_id_uint)?;
	}
	Ok(())
}

fn send_files(dir: &str, stream :&mut impl Write) -> Result<()>
{
	let dir = Path::new(dir);
	let n_files = read_dir(dir).unwrap().count() as u32;
    stream.write_u32::<NetworkEndian>(n_files)?;

    if dir.is_dir() {
        println!("Sending files...");

//		let mut rng = rand::thread_rng();
		let mut file_names = read_dir(dir)?.collect::<Vec<_>>();
//		file_names.as_mut_slice().shuffle(&mut rng);

//        for (_i, entry) in read_dir(dir)?.enumerate() {
		for (_i, entry) in file_names.into_iter().enumerate()
		{  
//            println!("File {}...", _i);
            let path = entry?.path();
            let mut f = File::open(path)?;

            let n_elems = f.seek(SeekFrom::End(0))? as usize/size_of::<u32>()-2;
            f.seek(SeekFrom::Start(0))?;

            let mut f = BufReader::new(f);
            send_single_file(&mut f, stream, n_elems)?;
        }
    } else {
        panic!("Empty directory.");
    }
    Ok(())
}

fn send_query_file(path: &str, stream :&mut impl Write) -> Result<()>
{
	println!("Sending query file ...");
	let path = Path::new(path);
	let mut f = File::open(path)?;

	let n_elems = f.seek(SeekFrom::End(0))? as usize / size_of::<u32>();
	f.seek(SeekFrom::Start(0))?;
	let mut f = BufReader::new(f);
	stream.write_u32::<NetworkEndian>(n_elems as u32)?;

	for _ in 0..n_elems
	{
		let rs_id_uint = f.read_u32::<NativeEndian>()?;
		stream.write_u32::<NetworkEndian>(rs_id_uint)?;
	}
	Ok(())
}

fn launch_enclave(enclave_runner: &str, enclave: &str) -> Child {
    Command::new(enclave_runner)
        .arg(enclave)
        .spawn()
        .expect("Failed to spawn binary.")
}

fn launch_binary(bin: &str) -> Child {
    Command::new(bin)
        .spawn()
        .expect("Failed to spawn binary.")
}

fn main()
{
    let enclave_runner = env::args().nth(2);
    let enclave = env::args().nth(3);

    // lauch enclave
    let mut child = match &enclave_runner {
        Some(enclave_runner) => {
            match &enclave {
                Some(enclave) => launch_enclave(&enclave_runner[..], &enclave[..]),
                None => launch_binary(&enclave_runner[..]),
            }
        }
        None => panic!("Incorrect input."),
    };

    // start listening to enclave
    let host = "localhost:1234";
    let listener = TcpListener::bind(&host).unwrap();
    let stream = listener.accept().unwrap().0;
    let mut stream = EncryptedWriter::with_capacity(TCP_BUFFER_SIZE, stream, &DUMMY_KEY);

    // send files
    let dir = &env::args().nth(1).unwrap()[..];

	// First Pass: For updating MCSK
    send_files(dir, &mut stream).expect("Send files error.");
    stream.flush().unwrap();

	// Second Pass: For updating the main sketch
    send_files(dir, &mut stream).expect("Send files error.");
    stream.flush().unwrap();

	// Query Pass
	let qpath = "/home/ckockan/chr1_uniq/chr1_uniq.ckz0";
    send_query_file(qpath, &mut stream).expect("Send files error.");
    stream.flush().unwrap();

	// Third Pass: For updating the RHHT-PCC
    send_files(dir, &mut stream).expect("Send files error.");
    stream.flush().unwrap();

    let ecode = child.wait().expect("failed to wait on child");
    assert!(ecode.success());
}
