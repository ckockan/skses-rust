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

pub struct CskfMap
{
	s_thres: i32,
    seed: (u32, u32, u32, u32),
    map: Vec<f32>,
}

impl CskfMap
{
    pub fn new(rng: &mut ThreadRng, params: &Parameters) -> Self
	{
        Self
		{
			s_thres: 200,
            seed: (rng.gen::<u32>(), rng.gen::<u32>(), rng.gen::<u32>(), rng.gen::<u32>()),
            map: vec![0.; params.csk_width],
        }
    }

    #[inline]
    pub fn update_pos_cskf(&mut self, item: u32, count: f32)
	{
		let pos = cal_hash(item as u64, self.seed.0 as u64, self.seed.1 as u64) & (self.map.len() - 1) as u32;
		let hash2 = cal_hash(item as u64, self.seed.2 as u64, self.seed.3 as u64);

		let count_ =
		if hash2 & 0x1 == 0
		{
			 -1.0 * count
		}
		else
		{
			count
		};
			
        let counter = &mut self.map[pos as usize];
		*counter = *counter + count_;
    }

	#[inline]
	pub fn cskf_get_median_odd(&self, item: u64) -> f32
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
		return self.map[pos as usize] * (sign as f32);
	}
}
