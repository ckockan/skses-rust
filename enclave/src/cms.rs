//use std::hash::{BuildHasher, Hasher};
use rand::{rngs::ThreadRng, Rng};
//use byteorder::{NativeEndian, ByteOrder};
use std::ops::Range;
use crate::parameters::*;


#[inline]
pub fn cal_hash(x: u64, a: u64, b: u64) -> u32
{
    let mut result: u64 = a * x + b;
	result = (result & 0x7FFFFFFF) + (result >> 31);

    if result >= 0x7FFFFFFF
	{
        result = result  - (0x7FFFFFFF as u64);
    }
    
	return result as u32;
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

pub struct CmsMap
{
	st_len: i64,
    seed: (u32, u32),
    map: &'static mut [i16; WIDTH],
}

impl CmsMap
{
    pub fn new(rng: &mut ThreadRng, map: &'static mut [i16; WIDTH]) -> Self
	{
        Self
		{
			st_len: 0,
            seed: (rng.gen::<u32>(), rng.gen::<u32>()),
            map
        }
    }


    #[inline]
    pub fn cal_pos(&self, item: u32) -> u32
	{
		cal_hash(item as u64, self.seed.0 as u64, self.seed.1 as u64) & (WIDTH - 1) as u32 
    }

    #[inline]
    pub fn update_pos(&mut self, pos: u32, count: i16)
	{
		self.st_len = self.st_len + (count as i64);

        let counter = &mut self.map[pos as usize];
        match (*counter).checked_add(count)
		{
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

	/*
	#[inline]
	pub fn cms_query_median_even(&self, item: u64) -> i16
	{
		let mut values: Vec<i16> = Vec::new();
		let mut median: i16;
		let mut hash: u32;
		let mut pos: u32;
		let mut log_width: u16 = __builtin_ctz(m_cms->width);

		for i in 0..m_cms->depth
		{
			hash = cal_hash(item, self.seed[i] as u64, self.seed[i + 1] as u64]);
			pos = hash & (WIDTH - 1);
			values.push(self.map[pos as usize]);

			// Guarantee unbiased query
			values[i] -= ((m_cms->st_length - values[i]) >> log_width);
		}

		// Sort values
		qsort(values, m_cms->depth, sizeof(int16_t), cmpfunc_int16);

		// Get median of values
		if(values[m_cms->depth / 2] < -(m_cms->s_thres))
		{
			median = values[m_cms->depth / 2 - 1];
		}
		else if(values[m_cms->depth / 2 - 1] > m_cms->s_thres)
		{
			median = values[m_cms->depth / 2];
		}
		else
		{
			median = (values[m_cms->depth / 2 - 1] + values[m_cms->depth / 2]) / 2;
		}
		return median;
	}
	*/
}
