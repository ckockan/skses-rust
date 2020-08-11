extern crate ndarray;

use ndarray::prelude::*;

pub fn matrix_ortho_proj(omat: ArrayView<f32, Ix2>, vec: ArrayView<f32, Ix1>, mut result: ArrayViewMut<f32, Ix1>, k: usize, m: usize)
{
	// In matrix form: result = omat^T * omat * vec
	for i in 0..m
	{
		for j in 0..m
		{
			for l in 0..k
			{
				result[i] += omat[[l, j]] * omat[[l, i]] * vec[j];
			}
		}
	}
}

pub fn matrix_ortho_proj1(mut omat: ArrayViewMut<f32, Ix2>, vec: ArrayView<f32, Ix1>, k: usize, m: usize)
{
	// In matrix form: result = omat^T * omat * vec
	for i in 0..m
	{
		for j in 0..m
		{
			for l in 0..k
			{
				omat[[k, i]] += omat[[l, j]] * omat[[l, i]] * vec[j];
			}
		}
	}
}
