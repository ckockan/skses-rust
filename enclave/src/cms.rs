//use std::hash::{BuildHasher, Hasher};
use rand::{rngs::ThreadRng, Rng};
//use byteorder::{NativeEndian, ByteOrder};
use std::ops::Range;
use crate::parameters::*;


#[inline]
pub fn cal_hash(x: u32, a: u32, b: u32) -> u32 {
    let mut result = a*x + b;
    if result >= 0x7FFFFFFF {
        result -= 0x7FFFFFFF;
    }
    result
}

//pub struct CmsHasher {
    //seed: (u64, u64),
    //current_result: u64,
//}

//impl CmsHasher {
    //pub fn new(seed: (u64, u64)) -> Self {
        //Self { seed, current_result: 0 }
    //}
//}

//impl Hasher for CmsHasher {
    //fn finish(&self) -> u64 {
        //self.current_result
    //}

    //fn write(&mut self, bytes: &[u8]) {
        //let x = NativeEndian::read_u32(&bytes[..4]);
        //self.current_result = cal_hash(x, self.seed.0, self.seed.1) as u64;
    //}
//}


//pub struct CmsBuildHasher {
    //seed: (u64, u64),
//}

//impl CmsBuildHasher {
    //pub fn new() -> Self {
        //let mut rng = thread_rng();
        //Self { seed: (rng.gen::<u64>(), rng.gen::<u64>()) }
    //}
//}

//impl BuildHasher for CmsBuildHasher {
    //type Hasher = CmsHasher;
    //fn build_hasher(&self) -> Self::Hasher {
        //CmsHasher::new(self.seed)
    //}
//}

pub struct CmsMap {
    seed: (u32, u32),
    map: &'static mut [i16; WIDTH],
}

impl CmsMap {
    pub fn new(rng: &mut ThreadRng, map: &'static mut [i16; WIDTH]) -> Self {
        Self {
            seed: (rng.gen::<u32>(), rng.gen::<u32>()),
            map
        }
    }


    #[inline]
    pub fn cal_pos(&self, item: u32) -> u32 {
        cal_hash(item, self.seed.0, self.seed.1) & 
            (WIDTH-1) as u32 
    }

    #[inline]
    pub fn update_pos(&mut self, pos: u32, count: i16) {
        let counter = &mut self.map[pos as usize];
        match (*counter).checked_add(count) {
            Some(c) => *counter = c,
            None => (),
        }
    }

    #[inline]
    pub fn update_pos_in_range(&mut self, pos: u32, count: i16, 
                               range: &Range<u32>) {
        if !range.contains(&pos) {
            return;
        }
        let counter = &mut self.map[pos as usize];
        match (*counter).checked_add(count) {
            Some(c) => *counter = c,
            None => (),
        }
    }
}
