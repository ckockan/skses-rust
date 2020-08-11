#[cfg(feature = "small-table")]
pub const WIDTH: usize = 0x40000;  // 0.25 MB

#[cfg(not(feature = "small-table"))]
pub const WIDTH: usize = 0x200000; // 2 MB

pub const DEPTH: usize = 8;

pub const MCSK_WIDTH: usize = 2000;
pub const MCSK_DEPTH: usize = 4096;
pub const MCSK_NUM_PC: usize = 2;

pub const N_THREAD: usize = 3;
pub const BUF_SIZE: usize = 1000; 
pub const N_FILE_QUOTA: usize = 6;
pub const PARTITION_SIZE: usize = 0x40000; 
pub const N_PARTITIONS: usize = (WIDTH+PARTITION_SIZE-1)/PARTITION_SIZE;
pub const TCP_BUFFER_SIZE: usize = 0x20000;
pub const DUMMY_KEY: [u8; 16] = [0u8; 16];
