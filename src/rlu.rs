use anyhow::{format_err, Result};
use std::iter::zip;

use crate::debug::debug;
use crate::dfs::DFS;
use crate::traits::{Int, Scalar};

// Simplified Compressed Column Storage
//
// > We shall represent a column vector as a sequence of records, each containing a
// value and a row index. The row indices need not be in increasing order. We shall
// represent a matrix as an array of column vectors indexed from 0 to n.
pub type Record<I, S> = (I, S);
pub type Col<I, S> = Vec<Record<I, S>>;
pub type Matrix<I, S> = Vec<Col<I, S>>;

pub fn solve<I: Int, S: Scalar, P: Int>(
    n: usize,
    a_rowidx: &[I],
    a_colptr: &[I],
    a_values: &[S],
    col_perm: Option<&[P]>,
    rhs: &mut [S],
    trans: bool,
) -> Result<()> {
    if rhs.len() % n != 0 {
        return Err(format_err!(
            "len rhs ({}) must be a multiple of n ({})",
            rhs.len(),
            n
        ));
    }

    let (l_mat, u_mat, p) =
        lu_decomposition(n, a_rowidx, a_colptr, a_values, col_perm, None, None, true);

    let mut x = vec![S::zero(); n];

    rhs.chunks_exact_mut(n).try_for_each(|b| {
        for i in 0..n {
            x[p[i].unwrap()] = b[i];
        }

        if !trans {
            lsolve(&l_mat, &mut x);
            usolve(&u_mat, &mut x);
        } else {
            utsolve(&u_mat, &mut x);
            ltsolve(&l_mat, &mut x);
        }

        match col_perm {
            Some(cperm) => {
                for i in 0..n {
                    b[cperm[i].to_index()] = x[i];
                }
            }
            None => b.copy_from_slice(&x),
        }

        Ok(())
    })
}

// 1. for j:= to n do
// 2.   {Compute column j of U and L.}
// 3.   Solve Ljuj = aj for uj;
// 4.   b'j := a'j - L'juj;
// 5.   Pivot: Swap bjj with the largest-magnitude element of b'j;
// 6.   ujj := bjj;
// 7.   l'j := b'j / ujj;
// 8. od

// LU decomposition (Gilbert-Peierls)
//
// Note: A is LU-decomposable <=> all principal minors are nonsingular
pub fn lu_decomposition<I: Int, S: Scalar, P: Int>(
    n: usize,
    a_rowidx: &[I],
    a_colptr: &[I], // n+1
    a_values: &[S],
    col_perm: Option<&[P]>,
    mut reachable: Option<&mut Vec<Vec<usize>>>,
    mut row_perm_inv: Option<&mut Vec<usize>>,
    pivot: bool,
) -> (Matrix<I, S>, Matrix<I, S>, Vec<Option<usize>>) {
    if let Some(reachable) = &reachable {
        assert_eq!(reachable.len(), 0);
    }
    if let Some(row_perm_inv) = &row_perm_inv {
        assert_eq!(row_perm_inv.len(), n);
    }

    let mut dfs = DFS::new(n);

    // row_perm(r) = Some(s) means row r of A is row s < jcol of PA (LU = PA).
    // row_perm(r) = None means row r of A has not yet been used as a
    // pivot and is therefore still below the diagonal.
    let mut row_perm: Vec<Option<usize>> = vec![None; n];

    let mut l_mat: Matrix<I, S> = vec![vec![]; n];
    let mut u_mat: Matrix<I, S> = vec![vec![]; n];

    // > We compute uj as a dense n-vector, so that in step 3.3 we can subtract a multiple
    // of column k of Lj from it in constant time per nonzero in that column.
    let mut x = vec![S::zero(); n];

    for k in 0..n {
        let kp = match col_perm {
            Some(perm) => perm[k].to_index(),
            None => k,
        };
        debug!("\nk = {}, kp = {}", k, kp);

        #[cfg(feature = "debug")]
        {
            print!("U =\n{}", crate::matrix_table(&u_mat));
            print!("L =\n{}", crate::matrix_table(&l_mat));
        }
        debug!("rperm = {:?}", row_perm);

        let b_rowidx = &a_rowidx[a_colptr[kp].to_index()..a_colptr[kp + 1].to_index()];
        let b_values = &a_values[a_colptr[kp].to_index()..a_colptr[kp + 1].to_index()];

        // Depth-first search from each nonzero of column jcol of A.
        let found = dfs.ludfs(&l_mat, b_rowidx, &row_perm);
        if let Some(reaches) = &mut reachable {
            reaches.push(found.to_vec());
        }

        // Compute the values of column jcol of L and U in the *dense* vector,
        // allocating storage for fill in L as necessary.
        lucomp(&l_mat, b_rowidx, b_values, &mut x, &row_perm, &found); // s = L\A[:,j]
        debug!("x = {:?}", x);

        let d = x[kp]; // TODO: numerically zero diagonal element at column check

        // Partial pivoting, diagonal elt. has max. magnitude in L.
        let (pivrow, maxabs) = found
            .iter()
            .filter(|i| row_perm[**i].is_none())
            .map(|&i| (i, x[i].norm()))
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .expect("pivot must exist");

        debug!("pivrow = {}, maxpiv = {:.6}, d = {:.6}", pivrow, maxabs, d);

        // TODO: Threshold pivoting.
        let pivt = if !pivot || (row_perm[kp].is_none() && d.norm() >= maxabs) {
            // No pivoting, diagonal element has irow = jcol.
            (kp, d)
        } else {
            (pivrow, x[pivrow])
        };

        // Copy the column elements of U, throwing out zeros.
        u_mat[k] = found
            .iter()
            .filter(|i| row_perm[**i].is_some())
            .map(|i| (I::from_usize(row_perm[*i].unwrap()), x[*i]))
            .collect();

        // Swap the max. value in L with the diagonal U(k,k).
        u_mat[k].push((I::from_usize(k), pivt.1));
        u_mat[k].shrink_to_fit(); // free up unused memory

        // Record the pivot in P.
        row_perm[pivt.0] = Some(k);

        if let Some(row_perm_inv) = &mut row_perm_inv {
            row_perm_inv[k] = pivt.0;
        }

        // Copy the column elements of L, throwing out zeros.
        l_mat[k] = found
            .iter()
            .filter(|i| row_perm[**i].is_none())
            .map(|i| (I::from_usize(*i), x[*i] / pivt.1))
            .collect();
        l_mat[k].shrink_to_fit(); // free up unused memory

        // > Since we know the nonzero structure of uj before we start,
        // we need only initialize and manipulate the positions in this
        // dense vector that correspond to nonzero positions.
        found.iter().for_each(|i| x[*i] = S::zero());

        // debug!("U[:,k] = {:?}", u_mat[k]);
        // debug!("L[:,k] = {:?}", l_mat[k]);
    }

    // Renumber the rows so the data structure represents L and U, not PtL and PtU.
    for row in &mut l_mat {
        for e in row {
            match row_perm[e.0.to_index()] {
                Some(e0) => e.0 = I::from_usize(e0),
                None => panic!(),
            }
        }
    }
    debug!("L =\n{}", crate::matrix_table(&l_mat));

    (l_mat, u_mat, row_perm)
}

// LU re-decomposition (Gilbert-Peierls).
//
// The pattern of the input matrix (`a_rowidx`, `a_colptr`) must be the same as the pattern
// passed to `lu_decomposition`. If `row_perm_inv` is specified pivoting will not be performed.
pub fn lu_redecomposition<I: Int, S: Scalar, P: Int>(
    n: usize,
    a_rowidx: &[I],
    a_colptr: &[I], // n+1
    a_values: &[S],
    reachable: &[Vec<usize>],
    row_perm_inv: Option<&[usize]>,
    col_perm: Option<&[P]>,
    pivot: bool,
) -> (Matrix<I, S>, Matrix<I, S>) {
    let mut row_perm: Vec<Option<usize>> = vec![None; n];

    let mut l_mat: Matrix<I, S> = vec![vec![]; n];
    let mut u_mat: Matrix<I, S> = vec![vec![]; n];

    // > We compute uj as a dense n-vector, so that in step 3.3 we can subtract a multiple
    // of column k of Lj from it in constant time per nonzero in that column.
    let mut x = vec![S::zero(); n];

    // let mut nz = 0;
    for k in 0..n {
        let kp = match col_perm {
            Some(perm) => perm[k].to_index(),
            None => k,
        };
        debug!("\nk = {}, kp = {}", k, kp);

        #[cfg(feature = "debug")]
        {
            print!("U =\n{}", crate::matrix_table(&u_mat));
            print!("L =\n{}", crate::matrix_table(&l_mat));
        }
        debug!("rperm = {:?}", row_perm);

        let b_rowidx = &a_rowidx[a_colptr[kp].to_index()..a_colptr[kp + 1].to_index()];
        let b_values = &a_values[a_colptr[kp].to_index()..a_colptr[kp + 1].to_index()];

        // Depth-first search from each nonzero of column jcol of A.
        let found = &reachable[k];
        debug!("found = {:?}", found);

        // Compute the values of column jcol of L and U in the *dense* vector,
        // allocating storage for fill in L as necessary.
        lucomp(&l_mat, b_rowidx, b_values, &mut x, &row_perm, found); // s = L\A[:,j]
        debug!("x = {:?}", x);

        let d = x[kp]; // TODO: numerically zero diagonal element at column check

        let pivt = match row_perm_inv {
            Some(row_perm_inv) => {
                let pivt0 = row_perm_inv[k];
                let pivt1 = if pivt0 == kp { d } else { x[pivt0] };
                (pivt0, pivt1)
            }
            None => {
                // Partial pivoting, diagonal elt. has max. magnitude in L.
                let (pivrow, maxabs) = found
                    .iter()
                    .filter(|i| row_perm[**i].is_none())
                    .map(|&i| (i, x[i].norm()))
                    .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                    .expect("pivot must exist");

                debug!("pivrow = {}, maxpiv = {:.6}, d = {:.6}", pivrow, maxabs, d);

                // TODO: Threshold pivoting.
                if !pivot || (row_perm[kp].is_none() && d.norm() >= maxabs) {
                    // No pivoting, diagonal element has irow = jcol.
                    (kp, d)
                } else {
                    (pivrow, x[pivrow])
                }
            }
        };

        // Copy the column elements of U, throwing out zeros.
        u_mat[k] = found
            .iter()
            .filter(|i| row_perm[**i].is_some())
            .map(|i| (I::from_usize(row_perm[*i].unwrap()), x[*i]))
            .collect();

        // Swap the max. value in L with the diagonal U(k,k).
        u_mat[k].push((I::from_usize(k), pivt.1));
        u_mat[k].shrink_to_fit(); // free up unused memory

        // Record the pivot in P.
        row_perm[pivt.0] = Some(k);

        // Copy the column elements of L, throwing out zeros.
        l_mat[k] = found
            .iter()
            .filter(|i| row_perm[**i].is_none())
            .map(|i| (I::from_usize(*i), x[*i] / pivt.1))
            .collect();
        l_mat[k].shrink_to_fit(); // free up unused memory

        // > Since we know the nonzero structure of uj before we start,
        // we need only initialize and manipulate the positions in this
        // dense vector that correspond to nonzero positions.
        found.iter().for_each(|i| x[*i] = S::zero());

        // debug!("U[:,k] = {:?}", u_mat[k]);
        // debug!("L[:,k] = {:?}", l_mat[k]);
    }

    // Renumber the rows so the data structure represents L and U, not PtL and PtU.
    for row in &mut l_mat {
        for e in row {
            match row_perm[e.0.to_index()] {
                Some(e0) => e.0 = I::from_usize(e0),
                None => panic!(),
            }
        }
    }
    debug!("L =\n{}", crate::matrix_table(&l_mat));

    (l_mat, u_mat)
}

fn lucomp<I: Int, S: Scalar>(
    l_mat: &Matrix<I, S>,
    b_rowidx: &[I],
    b_values: &[S],
    x: &mut Vec<S>,
    rperm: &[Option<usize>],
    found: &[usize],
) {
    // FOR(e, b) x[e->fst] = e->snd;
    // FOR(e, x) FOR(l, L[e->fst])
    //   x[l->fst] -= l->snd * e->snd;

    for (bi, bx) in zip(b_rowidx, b_values) {
        x[bi.to_index()] = *bx; // scatter
    }
    debug!("x = {:?}", x);

    // for e0 in 0..x.len() {
    for j in found {
        let e0 = match rperm[*j] {
            Some(jp) => jp,
            None => continue,
        };
        for l in &l_mat[e0] {
            let e1 = x[*j];

            x[l.0.to_index()] -= l.1 * e1;
        }
    }
}

pub fn lsolve<I: Int, S: Scalar>(l_mat: &Matrix<I, S>, b: &mut [S]) {
    // FOR(e, b) x[e->fst] = e->snd;
    // FOR(e, x) FOR(l, L[e->fst])
    //   x[l->fst] -= l->snd * e->snd;

    for e0 in 0..b.len() {
        for l in &l_mat[e0] {
            b[l.0.to_index()] -= l.1 * b[e0];
        }
    }
}

pub fn ltsolve<I: Int, S: Scalar>(l_mat: &Matrix<I, S>, b: &mut [S]) {
    for e0 in (0..b.len()).rev() {
        for l in l_mat[e0].iter().rev() {
            b[e0] -= l.1 * b[l.0.to_index()];
        }
    }
}

pub fn usolve<I: Int, S: Scalar>(u_mat: &Matrix<I, S>, b: &mut [S]) {
    // FORR(e, b) FORR(u, U[e->fst])
    //   if (u->fst == e->fst) e->snd /= u->snd;
    //   else b[u->fst] -= u->snd * e->snd;

    for e0 in (0..b.len()).rev() {
        for u in u_mat[e0].iter().rev() {
            if u.0.to_index() == e0 {
                b[e0] /= u.1;
            } else {
                b[u.0.to_index()] -= u.1 * b[e0];
            }
        }
    }
}

pub fn utsolve<I: Int, S: Scalar>(u_mat: &Matrix<I, S>, b: &mut [S]) {
    for e0 in 0..b.len() {
        for u in u_mat[e0].iter() {
            if u.0.to_index() == e0 {
                b[e0] /= u.1;
            } else {
                b[e0] -= u.1 * b[u.0.to_index()];
            }
        }
    }
}
