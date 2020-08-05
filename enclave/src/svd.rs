#![allow(clippy::many_single_char_names)]
use ndarray::{Array1, ArrayViewMut1, ArrayViewMut2};

/// a: mxn, w: n, v: nxn
pub fn svdcomp(mut a: ArrayViewMut2<f32>, mut w: ArrayViewMut1<f32>, mut v: ArrayViewMut2<f32>) {
    assert!(a.nrows() > a.ncols());
    assert!(w.len() == a.ncols());
    assert!(v.nrows() == v.ncols());
    assert!(v.nrows() == a.ncols());

    let m = a.nrows();
    let n = a.ncols();

    let mut i: usize;
    let mut l = 0usize;
    let mut nm = 0usize;
    let mut flag: usize;

    let mut c: f64;
    let mut f: f64;
    let mut h: f64;
    let mut s: f64;
    let mut x: f64;
    let mut y: f64;
    let mut z: f64;

    let mut anorm = 0f64;
    let mut g = 0f64;
    let mut scale = 0f64;
    let mut rv1 = Array1::<f64>::zeros(n);

    /* Householder reduction to bidiagonal form */
    for i in 0..n {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = 0.;
        scale = 0.;
        s = 0.;
        if i < m {
            for k in i..m {
                scale += (a[[k, i]] as f64).abs();
            }
            if scale != 0. {
                for k in i..m {
                    a[[k, i]] = (a[[k, i]] as f64 / scale) as f32;
                    s += a[[k, i]] as f64 * a[[k, i]] as f64;
                }
            }
            f = a[[i, i]] as f64;
            g = -sign(s.sqrt(), f);
            h = f * g - s;
            a[[i, i]] = (f - g) as f32;
            if i != n - 1 {
                for j in l..n {
                    s = 0.;
                    for k in i..m {
                        s += a[[k, i]] as f64 * a[[k, j]] as f64;
                    }
                    f = s / h;
                    for k in i..m {
                        a[[k, j]] += (f * a[[k, i]] as f64) as f32;
                    }
                }
            }
            for k in i..m {
                a[[k, i]] = (a[[k, i]] as f64 * scale) as f32;
            }
        }
        w[i] = (scale * g) as f32;

        /* right-hand reduction */
        g = 0.;
        scale = 0.;
        s = 0.;
        if i < m && i != n - 1 {
            for k in l..n {
                scale += (a[[i, k]] as f64).abs();
            }
            if scale != 0. {
                for k in l..n {
                    a[[i, k]] = (a[[i, k]] as f64 / scale) as f32;
                    s += a[[i, k]] as f64 * a[[i, k]] as f64;
                }
                f = a[[i, l]] as f64;
                g = -sign(s.sqrt(), f);
                h = f * g - s;
                a[[i, l]] = (f - g) as f32;
                for k in l..n {
                    rv1[k] = a[[i, k]] as f64 / h;
                }
                if i != m - 1 {
                    for j in l..m {
                        s = 0.;
                        for k in l..n {
                            s += a[[j, k]] as f64 * a[[i, k]] as f64;
                        }
                        for k in l..n {
                            a[[j, k]] += (s * rv1[k]) as f32;
                        }
                    }
                }
                for k in l..n {
                    a[[i, k]] = (a[[i, k]] as f64 * scale) as f32;
                }
            }
        }
        anorm = f64::max(anorm, (w[i] as f64).abs() + rv1[i].abs());
    }

    /* accumulate the right-hand transformation */
    for i in (0..n).rev() {
        if i < n - 1 {
            if g != 0. {
                /* double division to avoid underflow */
                for j in l..n {
                    v[[j, i]] = ((a[[i, j]] as f64 / a[[i, l]] as f64) / g) as f32;
                }
                for j in l..n {
                    s = 0.;
                    for k in l..n {
                        s += a[[i, k]] as f64 * v[[k, j]] as f64;
                    }
                    for k in l..n {
                        v[[k, j]] += (s * v[[k, i]] as f64) as f32;
                    }
                }
            }
            for j in l..n {
                v[[i, j]] = 0.;
                v[[j, i]] = 0.;
            }
        }
        v[[i, i]] = 1.;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for i in (0..n).rev() {
        l = i + 1;
        g = w[i] as f64;
        if i < n - 1 {
            for j in l..n {
                a[[i, j]] = 0.;
            }
        }
        if g != 0. {
            g = 1. / g;
            if i != n - 1 {
                for j in l..n {
                    s = 0.;
                    for k in l..m {
                        s += a[[k, i]] as f64 * a[[k, j]] as f64;
                    }
                    f = s / a[[i, i]] as f64 * g;
                    for k in i..m {
                        a[[k, j]] += (f * a[[k, i]] as f64) as f32;
                    }
                }
            }
            for j in i..m {
                a[[j, i]] = (a[[j, i]] as f64 * g) as f32;
            }
        } else {
            for j in i..m {
                a[[j, i]] = 0.;
            }
        }
        a[[i, i]] += 1.;
    }

    /* diagonalize the bidiagonal form */
    for k in (0..n).rev() {
        /* loop over singular values */
        let mut its = 0;
        loop {
            /* loop over allowed iterations */
            flag = 1;
            l = k;
            loop {
                /* test for splitting */
                if l != 0 {
                    nm = l - 1;
                }
                if (rv1[l].abs() + anorm - anorm).abs() < f64::EPSILON {
                    flag = 0;
                    break;
                }
                if l != 0 {
                    if ((w[nm] as f64).abs() + anorm - anorm).abs() < f64::EPSILON {
                        break;
                    }
                } else {
                    break;
                }
                if l==0 {
                    break;
                }
                l -= 1;
            }
            if flag != 0 {
                s = 1.;
                for i in l..(k + 1) {
                    f = s * rv1[i];
                    if (f.abs() + anorm - anorm).abs() > f64::EPSILON {
                        g = w[i] as f64;
                        h = pythag(f, g);
                        w[i] = h as f32;
                        h = 1. / h;
                        c = g * h;
                        s = -f * h;
                        for j in 0..m {
                            y = if l != 0 {
                                a[[j, nm]] as f64
                            } else {
                                0.
                            };
                            z = a[[j, i]] as f64;
                            if l != 0 {
                                a[[j, nm]] = (y * c + z * s) as f32;
                            } else {
                            }
                            a[[j, i]] = (z * c - y * s) as f32;
                        }
                    }
                }
            }
            z = w[k] as f64;
            if l == k {
                /* convergence */
                if z < 0. {
                    /* make singular value nonnegative */
                    w[k] = (-z) as f32;
                    for j in 0..n {
                        v[[j, k]] = -v[[j, k]];
                    }
                }
                break;
            }
            if its >= 30 {
                panic!("No convergence after 30 iterations")
            }

            /* shift from bottom 2 x 2 minor */
            x = w[l] as f64;
            if k != 0 {
                nm = k - 1;
                y = w[nm] as f64;
                g = rv1[nm];
            } else {
                y = 0.;
                g = 0.; 
            }
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2. * h * y);
            g = pythag(f, 1.);
            f = ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x;

            /* next QR transformation */
            c = 1.;
            s = 1.;
            for j in l..(nm + 1) {
                if k == 0 {
                    break;
                }
                i = j + 1;
                g = rv1[i];
                y = w[i] as f64;
                h = s * g;
                g *= c;
                z = pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for jj in 0..n {
                    x = v[[jj, j]] as f64;
                    z = v[[jj, i]] as f64;
                    v[[jj, j]] = (x * c + z * s) as f32;
                    v[[jj, i]] = (z * c - x * s) as f32;
                }
                z = pythag(f, h);
                w[j] = z as f32;
                if z != 0. {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for jj in 0..m {
                    y = a[[jj, j]] as f64;
                    z = a[[jj, i]] as f64;
                    a[[jj, j]] = (y * c + z * s) as f32;
                    a[[jj, i]] = (z * c - y * s) as f32;
                }
            }
            rv1[l] = 0.;
            rv1[k] = f;
            w[k] = x as f32;
            its += 1;
        }
    }
}

#[inline]
fn sign(a: f64, b: f64) -> f64 {
    if b >= 0. {
        a.abs()
    } else {
        -a.abs()
    }
}

fn pythag(a: f64, b: f64) -> f64 {
    let at = a.abs();
    let bt = b.abs();
    if at > bt {
        let ct = bt / at;
        at * (1. + ct * ct).sqrt()
    } else if bt > 0. {
        let ct = at / bt;
        bt * (1. + ct * ct).sqrt()
    } else {
        0.
    }
}


#[rustfmt::skip]
#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{arr1, arr2, Array1, Array2};

    #[test]
    fn svdcomp_test() {
        let mut a = arr2(&[
            [35.00, 46.00, 30.00, 76.00, 63.00],
            [5.00, 18.00, 98.00, 48.00, 12.00],
            [13.00, 0.00, 27.00, 11.00, 85.00],
            [44.00, 27.00, 30.00, 0.00, 89.00],
            [7.00, 51.00, 63.00, 87.00, 64.00],
            [33.00, 80.00, 24.00, 71.00, 75.00],
            [75.00, 6.00, 21.00, 5.00, 83.00],
        ]);

        let mut w = Array1::<f32>::zeros(a.ncols());
        let mut v = Array2::<f32>::zeros((a.ncols(), a.ncols()));
        svdcomp(a.view_mut(), w.view_mut(), v.view_mut());
        let a_ref = arr2(&[
            [-0.42200100, -0.11572283, -0.17700049, 0.49630406, 0.28182906],
            [-0.28050315, -0.50071543, -0.26249635, -0.17928235, -0.69868076],
            [-0.27553123, 0.30118147, 0.72184753, 0.22415136, -0.26279265],
            [-0.33965769, 0.42611924, 0.07237666, -0.57681906, -0.16821861],
            [-0.46736458, -0.38180637, 0.29115435, 0.18464789, 0.02459771],
            [-0.47677249, -0.08706811, -0.07983945, -0.45086846, 0.54377854],
            [-0.32589358, 0.55699855, -0.53133708, 0.31811792, -0.19687556],
        ]);
        let w_ref = arr1(&[270.92199707, 119.14959717, 42.28359604, 28.05463791, 77.31535339]);
        let v_ref = arr2(&[
            [-0.28844631, 0.43927738, -0.83686858, 0.15259080, -0.01417700, ],
            [-0.36012056, -0.21759605, -0.13336363, -0.73837823, 0.50987917, ],
            [-0.38944298, -0.38668016, -0.09708209, -0.22957274, -0.79792482, ],
            [-0.46030995, -0.55501896, -0.02616615, 0.61389172, 0.32018900, ],
            [-0.65081543, 0.54965430, 0.52130157, 0.04412302, -0.02484363, ]
        ]);
        assert!(a.iter().zip(a_ref.iter()).all(|(x, y)| (x-y).abs() < f32::EPSILON));
        assert!(w.iter().zip(w_ref.iter()).all(|(x, y)| (x-y).abs() < f32::EPSILON));
        assert!(v.iter().zip(v_ref.iter()).all(|(x, y)| (x-y).abs() < f32::EPSILON));
    }
}
