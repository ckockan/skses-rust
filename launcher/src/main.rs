use std::env;
use std::fs::{read_dir, File};
use std::path::Path;
use std::net::{TcpListener};
use std::process::{Command, Child};
use std::io::{Write, Read, BufWriter, BufReader, Result, SeekFrom, Seek};
use std::mem::size_of;
use byteorder::{ReadBytesExt, WriteBytesExt, NetworkEndian, NativeEndian};


fn send_single_file(file: &mut impl Read, stream: &mut impl Write, 
                    n_elems: usize) -> Result<()> {
    let patient_status = file.read_u32::<NativeEndian>()?;
    stream.write_u32::<NetworkEndian>(patient_status)?;
    let num_het_start = file.read_u32::<NativeEndian>()?;
    stream.write_u32::<NetworkEndian>(num_het_start)?;
    stream.write_u32::<NetworkEndian>(n_elems as u32)?;

    for _ in 0..n_elems {
        let rs_id_uint = file.read_u32::<NativeEndian>()?;
        stream.write_u32::<NetworkEndian>(rs_id_uint)?;
    }
    Ok(())
}

fn send_files(dir: &str, stream :&mut impl Write) -> Result<()> {
    let dir = Path::new(dir);
    let n_files = read_dir(dir).unwrap().count() as u32;
    stream.write_u32::<NetworkEndian>(n_files)?;

    if dir.is_dir() {
        println!("Sending files...");
        for (_i, entry) in read_dir(dir)?.enumerate() {
            //println!("File {}...", _i);
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

fn main() {
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
    let mut stream = BufWriter::with_capacity(0x100000, listener.accept().unwrap().0);

    // send files
    let dir = &env::args().nth(1).unwrap()[..];
    send_files(dir, &mut stream).expect("Send files error.");
    stream.flush().unwrap();
    let ecode = child.wait()
        .expect("failed to wait on child");
    assert!(ecode.success());
}
