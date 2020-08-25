use std::ops::Range;
use rand::{rngs::ThreadRng, Rng};
use ndarray::{Array2};
use crate::parameters::*;
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

pub struct Mcsk
{
	width: u32,
	depth: u32,
	k: u32,
	epsilon: f32,
    seed: (u32, u32, u32, u32),
    pub mcsk: Array2<f32>,
}

impl Mcsk
{
    pub fn new(rng: &mut ThreadRng, params: &Parameters) -> Self
	{
        Self
		{
			//m: MCSK_WIDTH as u32,
			width: params.mcsk_width as u32,
			depth: params.mcsk_depth as u32,
			k: 2,
			epsilon: (5.0 / params.mcsk_depth as f32).sqrt(),
            seed: (rng.gen::<u32>(), rng.gen::<u32>(), rng.gen::<u32>(), rng.gen::<u32>()),
            mcsk: Array2::zeros((params.mcsk_depth, params.mcsk_width + 1))
        }
    }


    #[inline]
    pub fn mcsk_update(&mut self, item: u32, col: usize, count: f32)
	{
		let mut hash: u32;
		let row: u32;
		let count_: f32;

		hash = cal_hash(item as u64, self.seed.0 as u64, self.seed.1 as u64);
		row = hash & (self.depth - 1) as u32;

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
		self.mcsk[[row as usize, self.width as usize]] = self.mcsk[[row as usize, self.width as usize]] + count_;
	}

	pub fn mcsk_mean_centering(&mut self)
	{
		let w: u32 = self.width;
		let d: u32 = self.depth;

		for r in 0..d
		{
			for c in 0..w
			{
				// Note: be aware of the potential precision issue with f32 vs. f64
				self.mcsk[[r as usize, c as usize]] -= self.mcsk[[r as usize, w as usize]] / self.width as f32;
			}
		}
	}

	pub fn mcsk_print(&self)
	{
		let w: u32 = self.width;
		let d: u32 = self.depth;

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
