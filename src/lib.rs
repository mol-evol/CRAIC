//! CRAIC acceleration core.
//!
//! A single, biology-agnostic primitive lives here: posterior "match"
//! probabilities from a 3-state pair-HMM (match / insert / delete with affine
//! gaps), computed by forward-backward in log space. Everything the workbench
//! needs about alignment uncertainty is derived from this matrix in Python:
//!   * the posterior explorer renders it directly,
//!   * column reliability = the mean posterior over a column's asserted pairs,
//!   * the built-in aligner does maximum-expected-accuracy traceback over it.
//!
//! Sequences arrive already integer-encoded (0..k-1) and the emission model
//! (joint k*k + background k) is supplied by the caller, so no substitution
//! matrix or alphabet is hard-coded in Rust.

use pyo3::prelude::*;

const NEG_INF: f64 = f64::NEG_INFINITY;

#[inline]
fn lse2(a: f64, b: f64) -> f64 {
    if a == NEG_INF {
        return b;
    }
    if b == NEG_INF {
        return a;
    }
    let m = if a > b { a } else { b };
    m + ((a - m).exp() + (b - m).exp()).ln()
}

#[inline]
fn lse3(a: f64, b: f64, c: f64) -> f64 {
    lse2(lse2(a, b), c)
}

/// Forward-backward posterior decoding of a pair-HMM.
///
/// Returns (flat posterior matrix of shape m*n, m, n) where element (i, j) is
/// P(residue i of `a` is homologous to residue j of `b`). The GIL is released
/// for the computation, so a long run does not freeze the interface.
#[pyfunction]
#[pyo3(signature = (a, b, k, joint, bg, delta=0.02, epsilon=0.5))]
fn pair_posteriors(
    py: Python<'_>,
    a: Vec<usize>,
    b: Vec<usize>,
    k: usize,
    joint: Vec<f64>,
    bg: Vec<f64>,
    delta: f64,
    epsilon: f64,
) -> PyResult<(Vec<f64>, usize, usize)> {
    // Defence in depth: craic.accel validates these too, but a bad delta here
    // would produce NaN log-transitions that propagate silently through every
    // posterior rather than failing.
    if !(delta > 0.0 && delta < 0.5) {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "delta must be in (0, 0.5); got {delta}"
        )));
    }
    if !(epsilon >= 0.0 && epsilon < 1.0) {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "epsilon must be in [0, 1); got {epsilon}"
        )));
    }

    // The dominant cost of every reliability run, every perturbation replicate
    // and the posterior explorer. gotoh_align and mea_align already release the
    // GIL; this one did not, so the most expensive kernel was the one that could
    // freeze the interface while it ran.
    Ok(py.allow_threads(move || {
        let m = a.len();
        let n = b.len();
        let w = n + 1;

        // log emissions
        let ljoint: Vec<f64> = joint.iter().map(|v| v.max(1e-300).ln()).collect();
        let lbg: Vec<f64> = bg.iter().map(|v| v.max(1e-300).ln()).collect();

        // log transitions
        let ld = delta.ln(); // M -> gap (open)
        let le = epsilon.ln(); // gap -> gap (extend)
        let lmm = (1.0 - 2.0 * delta).ln(); // M -> M
        let lgm = (1.0 - epsilon).ln(); // gap -> M

        let size = (m + 1) * w;
        let mut fm = vec![NEG_INF; size];
        let mut fi = vec![NEG_INF; size];
        let mut fd = vec![NEG_INF; size];
        fm[0] = 0.0;

        // ---- forward ----
        for i in 0..=m {
            for j in 0..=n {
                if i == 0 && j == 0 {
                    continue;
                }
                let idx = i * w + j;
                if i > 0 && j > 0 {
                    let e = ljoint[a[i - 1] * k + b[j - 1]];
                    let p = i.checked_sub(1).unwrap() * w + (j - 1);
                    fm[idx] = e + lse3(fm[p] + lmm, fi[p] + lgm, fd[p] + lgm);
                }
                if i > 0 {
                    let e = lbg[a[i - 1]];
                    let p = (i - 1) * w + j;
                    fi[idx] = e + lse2(fm[p] + ld, fi[p] + le);
                }
                if j > 0 {
                    let e = lbg[b[j - 1]];
                    let p = i * w + (j - 1);
                    fd[idx] = e + lse2(fm[p] + ld, fd[p] + le);
                }
            }
        }

        // partition function: any state may terminate
        let z = lse3(fm[m * w + n], fi[m * w + n], fd[m * w + n]);

        // ---- backward ----
        let mut bm = vec![NEG_INF; size];
        let mut bi = vec![NEG_INF; size];
        let mut bd = vec![NEG_INF; size];
        bm[m * w + n] = 0.0;
        bi[m * w + n] = 0.0;
        bd[m * w + n] = 0.0;

        for i in (0..=m).rev() {
            for j in (0..=n).rev() {
                if i == m && j == n {
                    continue;
                }
                let idx = i * w + j;
                let mut sm = NEG_INF;
                let mut si = NEG_INF;
                let mut sd = NEG_INF;
                if i < m && j < n {
                    let t = bm[(i + 1) * w + (j + 1)] + ljoint[a[i] * k + b[j]];
                    sm = lse2(sm, t + lmm);
                    si = lse2(si, t + lgm);
                    sd = lse2(sd, t + lgm);
                }
                if i < m {
                    let t = bi[(i + 1) * w + j] + lbg[a[i]];
                    sm = lse2(sm, t + ld);
                    si = lse2(si, t + le);
                }
                if j < n {
                    let t = bd[i * w + (j + 1)] + lbg[b[j]];
                    sm = lse2(sm, t + ld);
                    sd = lse2(sd, t + le);
                }
                bm[idx] = sm;
                bi[idx] = si;
                bd[idx] = sd;
            }
        }

        // ---- posterior match probabilities ----
        let mut post = vec![0.0f64; m * n];
        for i in 1..=m {
            for j in 1..=n {
                let v = fm[i * w + j] + bm[i * w + j] - z;
                let p = v.exp();
                post[(i - 1) * n + (j - 1)] = if p > 1.0 { 1.0 } else { p };
            }
        }
        (post, m, n)
    }))
}

/// Gotoh affine-gap alignment of two profiles, given their column score matrix.
///
/// `s` is the row-major `la * lb` match-score matrix (S[i, j] = score of aligning
/// column i of A with column j of B). Returns the aligned column indices for A and
/// B (an index, or -1 for a gap). This is the hot loop of the built-in aligner;
/// the GIL is released so the UI stays responsive while it runs.
#[pyfunction]
#[pyo3(signature = (s, la, lb, gap_open, gap_extend))]
fn gotoh_align(
    py: Python<'_>,
    s: Vec<f64>,
    la: usize,
    lb: usize,
    gap_open: f64,
    gap_extend: f64,
) -> PyResult<(Vec<i64>, Vec<i64>)> {
    Ok(py.allow_threads(move || {
        let neg = -1e18_f64;
        let w = lb + 1;
        let size = (la + 1) * w;
        let mut mm = vec![neg; size]; // match state
        let mut ix = vec![neg; size]; // gap in B (consume A)
        let mut iy = vec![neg; size]; // gap in A (consume B)
        let mut tbm = vec![0u8; size]; // traceback pointer for each state
        let mut tbx = vec![0u8; size];
        let mut tby = vec![0u8; size];
        mm[0] = 0.0;
        for i in 1..=la {
            ix[i * w] = gap_open + (i as f64 - 1.0) * gap_extend;
        }
        for j in 1..=lb {
            iy[j] = gap_open + (j as f64 - 1.0) * gap_extend;
        }
        for i in 1..=la {
            for j in 1..=lb {
                let sij = s[(i - 1) * lb + (j - 1)];
                let idx = i * w + j;
                // match: best of the three states on the diagonal, then add s
                let pd = (i - 1) * w + (j - 1);
                let (mut b, mut kb) = (mm[pd], 0u8);
                if ix[pd] > b { b = ix[pd]; kb = 1; }
                if iy[pd] > b { b = iy[pd]; kb = 2; }
                mm[idx] = b + sij;
                tbm[idx] = kb;
                // gap in B: from the cell above (i-1, j)
                let pu = (i - 1) * w + j;
                let (c0, c1, c2) = (mm[pu] + gap_open, ix[pu] + gap_extend, iy[pu] + gap_open);
                let (mut b, mut kb) = (c0, 0u8);
                if c1 > b { b = c1; kb = 1; }
                if c2 > b { b = c2; kb = 2; }
                ix[idx] = b;
                tbx[idx] = kb;
                // gap in A: from the cell to the left (i, j-1)
                let pl = i * w + (j - 1);
                let (c0, c1, c2) = (mm[pl] + gap_open, ix[pl] + gap_open, iy[pl] + gap_extend);
                let (mut b, mut kb) = (c0, 0u8);
                if c1 > b { b = c1; kb = 1; }
                if c2 > b { b = c2; kb = 2; }
                iy[idx] = b;
                tby[idx] = kb;
            }
        }
        let e = la * w + lb;
        let mut b = mm[e];
        let mut state = 0u8;
        if ix[e] > b { b = ix[e]; state = 1; }
        if iy[e] > b { state = 2; }
        let _ = b;
        let mut i = la;
        let mut j = lb;
        let mut cols_a: Vec<i64> = Vec::with_capacity(la + lb);
        let mut cols_b: Vec<i64> = Vec::with_capacity(la + lb);
        while i > 0 || j > 0 {
            let idx = i * w + j;
            if i > 0 && j > 0 && state == 0 {
                cols_a.push((i - 1) as i64);
                cols_b.push((j - 1) as i64);
                state = tbm[idx];
                i -= 1;
                j -= 1;
            } else if i > 0 && (state == 1 || j == 0) {
                cols_a.push((i - 1) as i64);
                cols_b.push(-1);
                state = tbx[idx];
                i -= 1;
            } else {
                cols_a.push(-1);
                cols_b.push((j - 1) as i64);
                state = tby[idx];
                j -= 1;
            }
        }
        cols_a.reverse();
        cols_b.reverse();
        (cols_a, cols_b)
    }))
}

/// Pure maximum-expected-accuracy (gap-free) alignment of a posterior score matrix.
///
/// `s` is the row-major `la * lb` posterior "match" score matrix. Unlike the Gotoh
/// aligner there is no gap penalty: the alignment maximises the expected number of
/// correctly-aligned residue pairs (ProbCons decoding). Returns the aligned column
/// indices for A and B (an index, or -1 for a gap). Mirrors `_mea_numpy` in
/// accel.py exactly (cross-validated in the test-suite). GIL released while it runs.
#[pyfunction]
#[pyo3(signature = (s, la, lb))]
fn mea_align(
    py: Python<'_>,
    s: Vec<f64>,
    la: usize,
    lb: usize,
) -> PyResult<(Vec<i64>, Vec<i64>)> {
    Ok(py.allow_threads(move || {
        let w = lb + 1;
        let size = (la + 1) * w;
        let mut mm = vec![0.0f64; size]; // expected-accuracy grid; borders stay 0
        let mut tb = vec![0u8; size]; // 0 diag, 1 up (gap in B), 2 left (gap in A)
        for i in 1..=la {
            for j in 1..=lb {
                let diag = mm[(i - 1) * w + (j - 1)] + s[(i - 1) * lb + (j - 1)];
                let up = mm[(i - 1) * w + j];
                let left = mm[i * w + (j - 1)];
                let idx = i * w + j;
                if diag >= up && diag >= left {
                    mm[idx] = diag;
                    tb[idx] = 0;
                } else if up >= left {
                    mm[idx] = up;
                    tb[idx] = 1;
                } else {
                    mm[idx] = left;
                    tb[idx] = 2;
                }
            }
        }
        let mut i = la;
        let mut j = lb;
        let mut cols_a: Vec<i64> = Vec::with_capacity(la + lb);
        let mut cols_b: Vec<i64> = Vec::with_capacity(la + lb);
        while i > 0 || j > 0 {
            let t = if i > 0 && j > 0 {
                tb[i * w + j]
            } else if j == 0 {
                1
            } else {
                2
            };
            if t == 0 {
                cols_a.push((i - 1) as i64);
                cols_b.push((j - 1) as i64);
                i -= 1;
                j -= 1;
            } else if t == 1 {
                cols_a.push((i - 1) as i64);
                cols_b.push(-1);
                i -= 1;
            } else {
                cols_a.push(-1);
                cols_b.push((j - 1) as i64);
                j -= 1;
            }
        }
        cols_a.reverse();
        cols_b.reverse();
        (cols_a, cols_b)
    }))
}

#[pymodule]
fn _accel(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pair_posteriors, m)?)?;
    m.add_function(wrap_pyfunction!(gotoh_align, m)?)?;
    m.add_function(wrap_pyfunction!(mea_align, m)?)?;
    m.add("__doc__", "CRAIC acceleration core (pair-HMM posteriors + Gotoh / MEA aligners).")?;
    Ok(())
}
