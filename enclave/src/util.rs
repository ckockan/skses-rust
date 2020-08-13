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

pub fn orthonormal_test(V: ArrayView<f32, Ix2>, m: usize)
{
	// Test whether the row vectors are orthogonal
	let mut dot12: f32 = 0.0;
	let mut dot13: f32 = 0.0;
	let mut dot23: f32 = 0.0;
	let mut norm1: f32 = 0.0;
	let mut norm2: f32 = 0.0;
	let mut norm3: f32 = 0.0;

	for i in 0..m
	{
		dot12 += V[[i, 0]] * V[[i, 1]];
		dot13 += V[[i, 0]] * V[[i, 2]];
		dot23 += V[[i, 1]] * V[[i, 2]];
		norm1 += V[[i, 0]] * V[[i, 0]];
		norm2 += V[[i, 1]] * V[[i, 1]];
		norm3 += V[[i, 2]] * V[[i, 2]];
	}
	println!("{}\t{}\t{}\n", dot12, dot13, dot23);
	println!("{}\t{}\t{}\n", norm1, norm2, norm3);
}

pub fn orthonormal_test_t(V: ArrayView<f32, Ix2>, m: usize)
{
	let mut dot12: f32 = 0.0;
	let mut dot13: f32 = 0.0;
	let mut dot23: f32 = 0.0;
	let mut norm1: f32 = 0.0;
	let mut norm2: f32 = 0.0;
	let mut norm3: f32 = 0.0;

	for i in 0..m
	{
		dot12 += V[[0, i]] * V[[1, i]];
		dot13 += V[[0, i]] * V[[2, i]];
		dot23 += V[[1, i]] * V[[2, i]];
		norm1 += V[[0, i]] * V[[0, i]];
		norm2 += V[[1, i]] * V[[1, i]];
		norm3 += V[[2, i]] * V[[2, i]];
	}
	println!("{}\t{}\t{}\n", dot12, dot13, dot23);
	println!("{}\t{}\t{}\n", norm1, norm2, norm3);
}
