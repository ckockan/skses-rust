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

pub struct CskMap
{
	s_thres: i32,
    seed: (u32, u32, u32, u32),
    map: &'static mut [i16; WIDTH],
}

impl CskMap
{
    pub fn new(rng: &mut ThreadRng, map: &'static mut [i16; WIDTH]) -> Self
	{
        Self
		{
			s_thres: 200,
            seed: (rng.gen::<u32>(), rng.gen::<u32>(), rng.gen::<u32>(), rng.gen::<u32>()),
            map
        }
    }

	/*
    #[inline]
    pub fn cal_pos(&self, item: u32) -> u32
	{
		let mut hash1 = cal_hash(item as u64, self.seed.0 as u64, self.seed.1 as u64) & (WIDTH - 1) as u32;
		let mut hash2 = cal_hash(item as u64, self.seed.2 as u64, self.seed.3 as u64);
		let count_ = ((hash2 & 0x1) == 0) > -1 : 1) * count;
		return count_;
    }
	*/

    #[inline]
    pub fn update_pos_csk(&mut self, item: u32, count: i16)
	{
		let mut pos = cal_hash(item as u64, self.seed.0 as u64, self.seed.1 as u64) & (WIDTH - 1) as u32;
		let mut hash2 = cal_hash(item as u64, self.seed.2 as u64, self.seed.3 as u64);

		let count_ =
		if hash2 & 0x1 == 0
		{
			 -1 * count
		}
		else
		{
			count
		};
			
        let counter = &mut self.map[pos as usize];
		/*
		if counter >= HASH_MAX_16 && count_ > 0
		{
			return;
		}
		if counter <= HASH_MIN_16 && count_ < 0
		{
			return;
		}
		*/

        match (*counter).checked_add(count_)
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

	#[inline]
	pub fn csk_get_median_odd(&self, item: u64) -> i16
	{
		let mut hash: u32;
		let mut pos: u32;
		let mut sign: i32;

		hash = cal_hash(item, self.seed.0 as u64, self.seed.1 as u64);
		pos = hash & (WIDTH - 1) as u32;
		hash = cal_hash(item, self.seed.2 as u64, self.seed.3 as u64);
		
		let sign =
		if hash & 0x1 == 0
		{
			-1
		}
		else
		{
			1
		};
		return self.map[pos as usize] * sign;
	}

	/*
	pub fn csk_query_median_odd(&self, item: u64) -> i16
	{
		let mut values = Vec::new();
		let mut median: i16;
		let mut hash: u32;
		let mut pos: u32;
		let mut sign: i32;

		hash = cal_hash(item, m_csk->seeds[i << 1], m_csk->seeds[(i << 1) + 1]);
			pos = hash & m_csk->width_minus_one;
			hash = cal_hash(item, m_csk->seeds[(i + m_csk->depth) << 1], m_csk->seeds[((i + m_csk->depth) << 1) + 1]);
			sign = ((hash & 0x1) == 0) ? -1 : 1;
			values[i] = m_csk->sketch[i][pos] * sign;

		// Sort values
		qsort(values, m_csk->depth, sizeof(int16_t), cmpfunc_int16);

		// Get median of values
		median = values[m_csk->depth / 2];

		// Free memory
		free(values);

		return median;
	}
	*/
}
