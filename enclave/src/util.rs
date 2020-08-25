use ndarray::{ArrayView1, ArrayView2, ArrayViewMut1, ArrayViewMut2};

pub fn matrix_ortho_proj(omat: ArrayView2<f32>, vec: ArrayView1<f32>, mut result: ArrayViewMut1<f32>, k: usize, m: usize)
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

pub fn matrix_ortho_proj1(mut omat: ArrayViewMut2<f32>, vec: ArrayView1<f32>, k: usize, m: usize)
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

pub fn orthonormal_test(v: ArrayView2<f32>, m: usize)
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
		dot12 += v[[i, 0]] * v[[i, 1]];
		dot13 += v[[i, 0]] * v[[i, 2]];
		dot23 += v[[i, 1]] * v[[i, 2]];
		norm1 += v[[i, 0]] * v[[i, 0]];
		norm2 += v[[i, 1]] * v[[i, 1]];
		norm3 += v[[i, 2]] * v[[i, 2]];
	}
	println!("{}\t{}\t{}\n", dot12, dot13, dot23);
	println!("{}\t{}\t{}\n", norm1, norm2, norm3);
}

pub fn orthonormal_test_t(v: ArrayView2<f32>, m: usize)
{
	let mut dot12: f32 = 0.0;
	let mut dot13: f32 = 0.0;
	let mut dot23: f32 = 0.0;
	let mut norm1: f32 = 0.0;
	let mut norm2: f32 = 0.0;
	let mut norm3: f32 = 0.0;

	for i in 0..m
	{
		dot12 += v[[0, i]] * v[[1, i]];
		dot13 += v[[0, i]] * v[[2, i]];
		dot23 += v[[1, i]] * v[[2, i]];
		norm1 += v[[0, i]] * v[[0, i]];
		norm2 += v[[1, i]] * v[[1, i]];
		norm3 += v[[2, i]] * v[[2, i]];
	}
	println!("{}\t{}\t{}\n", dot12, dot13, dot23);
	println!("{}\t{}\t{}\n", norm1, norm2, norm3);
}

/*
let fileQ = File::open("/home/ckockan/matrixV_norm_4096.out").unwrap();
		let readerQ = BufReader::new(fileQ);
		let mut qj = 0;
		for line in readerQ.lines()
		{
			let s = line.unwrap();
			//println!("{}", line.unwrap());
			let vecQ: Vec<&str> = s.split_whitespace().collect();
			for qi in 0..vecQ.len()
			{
//				println!("{}\t{}", qi, qj);
				Q[[qi, qj]] = vecQ[qi].parse().unwrap();
			}
			qj = qj + 1;
		}
//		println!("{:?}", Q);
*/
