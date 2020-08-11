#[inline]
pub fn cat_chi_sq(ssqg: u16, total: u16, dotprod: f32, sx: f32, sy: f32, sy2: f32, pc_projections: &Vec<f32>) -> f32
{
	let mut chi_sq_val: f32 = 0.0;
	let mut sx2: f32 = ssqg as f32;

	for pc in 0..2 //MCSK_NUM_PC
	{
		sx2 = sx2 - (pc_projections[pc] * pc_projections[pc]);
	}

	chi_sq_val = (total as f32 * dotprod) - (sx * sy);
	chi_sq_val = (chi_sq_val * chi_sq_val) * (total as f32) / (((total as f32) * sx2 - sx * sx) * ((total as f32) * sy2 - sy * sy));
	return chi_sq_val;
}

#[inline]
pub fn chi_sq(case_min: u16, control_min: u16, case_total: u16, control_total: u16) -> f32
{
	let case_maj: u16 = case_total - case_min;
	let control_maj: u16 = control_total - control_min;
	let pop_total: u16 = case_total + control_total;

	// Compute the observed frequencies
	let case_maj_f: f32 = (case_maj as f32) / (pop_total as f32);
	let case_min_f: f32 = (case_min as f32) / (pop_total as f32);
	let control_maj_f: f32 = (control_maj as f32) / (pop_total as f32);
	let control_min_f: f32 = (control_min as f32) / (pop_total as f32);
	let obs_freq = vec![case_maj_f, case_min_f, control_maj_f, control_min_f];

	// Compute the expected frequencies
	let case_total_f: f32 = (case_total as f32) / (pop_total as f32);
	let control_total_f: f32 = (control_total as f32) / (pop_total as f32);
	let maj_total_f: f32 = case_maj_f + control_maj_f;
	let min_total_f: f32 = case_min_f + control_min_f;
	let exp_freq = vec![case_total_f * maj_total_f,
					    case_total_f * min_total_f,
						control_total_f * maj_total_f,
						control_total_f * min_total_f];

	// Compute expected counts
	let exp_count = vec![exp_freq[0] * (pop_total as f32),
						 exp_freq[1] * (pop_total as f32),
						 exp_freq[2] * (pop_total as f32),
						 exp_freq[3] * (pop_total as f32)];

	// Compute the Chi-Squared value
	let mut chi_sq_val: f32 = 0.0;
	chi_sq_val = chi_sq_val + (((case_maj as f32) - exp_count[0]) * ((case_maj as f32) - exp_count[0]) / exp_count[0]);
	chi_sq_val = chi_sq_val + (((case_min as f32) - exp_count[1]) * ((case_min as f32) - exp_count[1]) / exp_count[1]);
	chi_sq_val = chi_sq_val + (((control_maj as f32) - exp_count[2]) * ((control_maj as f32) - exp_count[2]) / exp_count[2]);
	chi_sq_val = chi_sq_val + (((control_min as f32) - exp_count[3]) * ((control_min as f32) - exp_count[3]) / exp_count[3]);

	return chi_sq_val;	
}

pub fn poz(z: f32) -> f32
{
	let mut x: f32 = 0.0;
	if z != 0.0
	{
		let mut y: f32 = 0.5 * z.abs();
		if y >= 3.0
		{
			x = 1.0;
		}
		else if y < 1.0
		{
			let w: f32 = y * y;
			x = ((((((((0.000124818987 * w
				-0.001075204047) * w + 0.005198775019) * w
				-0.019198292004) * w + 0.059054035642) * w
				-0.151968751364) * w + 0.319152932694) * w
				-0.531923007300) * w + 0.797884560593) * y * 2.0;
		}
		else
		{
			y = y - 2.0;
			x = (((((((((((((-0.000045255659 * y
				+0.000152529290) * y -0.000019538132) * y
				-0.000676904986) * y +0.001390604284) * y
				-0.000794620820) * y -0.002034254874) * y
				+0.006549791214) * y -0.010557625006) * y
				+0.011630447319) * y -0.009279453341) * y
				+0.005353579108) * y -0.002141268741) * y
				+0.000535310849) * y +0.999936657524;
		}
	}

	if z > 0.0
	{
		return (x + 1.0) * 0.5;
	}
	else
	{
		return (1.0 - x) * 0.5;
	}
}

pub fn pochisq(x: f32) -> f32
{
	if x <= 0.0
	{
		return 1.0;
	}
	else
	{
		return 2.0 * poz(-x.sqrt());
	}
}
