extern crate ndarray;

use ndarray::prelude::*;
use rand::{rngs::ThreadRng, Rng};
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

pub struct Mcsk
{
	m: u32,
	k: u32,
	epsilon: f32,
    seed: (u32, u32, u32, u32),
    pub mcsk: Array<f32, Ix2>,
}

impl Mcsk
{
    pub fn new(rng: &mut ThreadRng) -> Self
	{
        Self
		{
			m: MCSK_WIDTH as u32,
			k: 2,
			epsilon: (5.0 / MCSK_DEPTH as f32).sqrt(),
            seed: (rng.gen::<u32>(), rng.gen::<u32>(), rng.gen::<u32>(), rng.gen::<u32>()),
            mcsk: Array::zeros((MCSK_DEPTH, MCSK_WIDTH + 1))
        }
    }


    #[inline]
    pub fn mcsk_update(&mut self, item: u32, col: usize, count: f32)
	{
		let mut hash: u32;
		let mut row: u32;
		let mut count_: f32;

		hash = cal_hash(item as u64, self.seed.0 as u64, self.seed.1 as u64);
		row = hash & (MCSK_DEPTH - 1) as u32;

		hash = cal_hash(item as u64, self.seed.2 as u64, self.seed.3 as u64);
		if hash & 0x1 == 0
		{
			count_ = -count;
		}
		else
		{
			count_ = count;
		}

		self.mcsk[[row as usize, col as usize]] = self.mcsk[[row as usize, col as usize]] + count_;
		self.mcsk[[row as usize, self.m as usize]] = self.mcsk[[row as usize, self.m as usize]] + count_;
	}

	pub fn mcsk_mean_centering(&mut self)
	{
		let mut d: u32 = MCSK_DEPTH as u32;
		let mut w: u32 = self.m;

		for r in 0..d
		{
			for c in 0..w
			{
				// Note: be aware of the potential precision issue with f32 vs. f64
				self.mcsk[[r as usize, c as usize]] -= self.mcsk[[r as usize, w as usize]] / self.m as f32;
			}
		}
	}

	pub fn mcsk_print(&self)
	{
		let mut d: u32 = MCSK_DEPTH as u32;
		let mut w: u32 = self.m;

		for r in 0..d
		{
			for c in 0..w
			{
				print!("{:.10} ", self.mcsk[[r as usize, c as usize]]);
			}
			println!();
			if r >= 10
			{
				break;
			}
		}
	}
}
