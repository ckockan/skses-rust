use std::hash::{BuildHasher, Hasher};
use rand::{rngs::ThreadRng, Rng, thread_rng};
use byteorder::{NativeEndian, ByteOrder};
use std::ops::Range;

#[inline]
pub fn cal_hash(x: u32, a: u64, b: u64) -> u64 {
    let mut result = a*x as u64 + b;
    if result >= 0x7FFFFFFF {
        result -= 0x7FFFFFFF;
    }
    result
}

pub struct CmsHasher {
    seed: (u64, u64),
    current_result: u64,
}

impl CmsHasher {
    pub fn new(seed: (u64, u64)) -> Self {
        Self { seed, current_result: 0 }
    }
}

impl Hasher for CmsHasher {
    fn finish(&self) -> u64 {
        self.current_result
    }

    fn write(&mut self, bytes: &[u8]) {
        let x = NativeEndian::read_u32(&bytes[..4]);
        self.current_result = cal_hash(x, self.seed.0, self.seed.1);
    }
}


pub struct CmsBuildHasher {
    seed: (u64, u64),
}

impl CmsBuildHasher {
    pub fn new() -> Self {
        let mut rng = thread_rng();
        Self { seed: (rng.gen::<u64>(), rng.gen::<u64>()) }
    }
}

impl BuildHasher for CmsBuildHasher {
    type Hasher = CmsHasher;
    fn build_hasher(&self) -> Self::Hasher {
        CmsHasher::new(self.seed)
    }
}

pub struct CmsMap {
    seed: (u64, u64),
    map: Vec<i16>,
    width_m_1: usize,
}

impl CmsMap {
    pub fn new(width: usize, rng: &mut ThreadRng) -> Self {
        Self {
            seed: (rng.gen::<u64>(), rng.gen::<u64>()),
            map: vec![0i16; width] ,
            width_m_1: width-1,
        }
    }

    #[inline]
    pub fn update(&mut self, item: u32, count: i16) {
        let pos = cal_hash(item, self.seed.0, self.seed.1) as usize & 
                   self.width_m_1;
        let counter = &mut self.map[pos];
        match (*counter).checked_add(count) {
            Some(c) => *counter = c,
            None => (),
        }
    }

    #[inline]
    pub fn update_in_range(&mut self, item: u32, count: i16, range: Range<usize>) {
        let pos = cal_hash(item, self.seed.0, self.seed.1) as usize & 
                   self.width_m_1;
        if !range.contains(&pos) {
            return;
        }
        let counter = &mut self.map[pos];
        match (*counter).checked_add(count) {
            Some(c) => *counter = c,
            None => (),
        }
    }
}
