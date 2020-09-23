use std::ops::Range;
use rand::{rngs::ThreadRng, Rng};
use crate::Parameters;

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
    map: Vec<i16>,
}

impl CskMap
{
    pub fn new(rng: &mut ThreadRng, params: &Parameters) -> Self
	{
        Self
		{
			s_thres: 200,
            seed: (rng.gen::<u32>(), rng.gen::<u32>(), rng.gen::<u32>(), rng.gen::<u32>()),
            map: vec![0; params.csk_width],
        }
    }

    #[inline]
    pub fn update_pos_csk(&mut self, item: u32, count: i16)
	{
		let pos = cal_hash(item as u64, self.seed.0 as u64, self.seed.1 as u64) & (self.map.len() - 1) as u32;
		let hash2 = cal_hash(item as u64, self.seed.2 as u64, self.seed.3 as u64);

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
        match (*counter).checked_add(count_)
		{
            Some(c) => *counter = c,
            None => (),
        }
    }

    #[inline]
    pub fn update_pos_in_range(&mut self, pos: u32, count: i16, range: &Range<u32>)
	{
        if !range.contains(&pos)
		{
            return;
        }

        let counter = &mut self.map[pos as usize];
        match (*counter).checked_add(count)
		{
            Some(c) => *counter = c,
            None => (),
        }
    }

	#[inline]
	pub fn csk_get_median_odd(&self, item: u64) -> i16
	{
		let mut hash: u32;
		let pos: u32;

		hash = cal_hash(item, self.seed.0 as u64, self.seed.1 as u64);
		pos = hash & (self.map.len() - 1) as u32;
		hash = cal_hash(item, self.seed.2 as u64, self.seed.3 as u64);
		
		let sign: i32 =
		if hash & 0x1 == 0
		{
			-1
		}
		else
		{
			1
		};
		return self.map[pos as usize] * (sign as i16);
	}
}
