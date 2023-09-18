use std::{
    cell::RefCell,
    collections::{BTreeMap, HashSet},
};

use sorted_vec::SortedVec;

// Simplified Compressed Column Storage
//
// > We shall represent a column vector as a sequence of records, each containing a
// value and a row index. The row indices need not be in increasing order. We shall
// represent a matrix as an array of column vectors indexed from 0 to n.
pub type Record = (usize, f64);
pub type Col = Vec<Record>;
pub type Matrix = Vec<Col>;

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
pub fn lu_decomposition(
    a_mat: &Matrix,
    col_perm: Option<&[usize]>,
    pivot: bool,
) -> (Matrix, Matrix, Vec<usize>) {
    let n = a_mat.len();

    // row_perm(r) = Some(s) means row r of A is row s < jcol of PA (LU = PA).
    // row_perm(r) = None means row r of A has not yet been used as a
    // pivot and is therefore still below the diagonal.
    let mut pivots = HashSet::new();
    // let mut row_perm: Vec<Option<usize>> = vec![None; n];
    // let mut row_perm_inv: Vec<Option<usize>> = vec![None; n];
    let mut p: Vec<usize> = (0..n).collect();

    let mut l_mat: Matrix = vec![vec![]; n];
    let mut u_mat: Matrix = vec![vec![]; n];

    // > We compute uj as a dense n-vector, so that in step 3.3 we can subtract a multiple
    // of column k of Lj from it in constant time per nonzero in that column.
    let mut x = vec![0.0; n];

    for k in 0..n {
        let kp = match col_perm {
            Some(perm) => perm[k],
            None => k,
        };
        // Compute the values of column jcol of L and U in the *dense* vector,
        // allocating storage for fill in L as necessary.
        // let mut s: BTreeMap<usize, RefCell<f64>> = sp_lsolve(&l_mat, &a_mat[kp], &p); // s = L\A[:,j]
        sp_lsolve_d(&l_mat, &a_mat[kp], &mut x); // s = L\A[:,j]
                                                 // assert_ne!(s.len(), 0);

        // let d = match s.get(&k) {
        //     Some(d) => d.borrow().to_owned(),
        //     // None => s.last_key_value().unwrap().1.borrow().to_owned(),
        //     None => 0.0, // numerically zero diagonal element at column
        // };
        let d = x[k];

        // No pivoting, diagonal element has irow = jcol.
        // Partial pivoting, diagonal elt. has max. magnitude in L.
        // let (pivrow, maxpiv) = s
        //     .range(k + 1..)
        //     .max_by(|(_, v0_rc), (_, v1_rc)| {
        //         let v0 = v0_rc.borrow().abs();
        //         let v1 = v1_rc.borrow().abs();
        //         v0.partial_cmp(&v1).unwrap_or(std::cmp::Ordering::Equal)
        //     })
        //     .map(|(i, v_rc)| (*i, v_rc.borrow().to_owned()))
        //     .unwrap_or((k, d));
        let mut pivrow = k;
        let mut maxpiv = d;
        for (i, xi) in x.iter().enumerate() {
            if i > k {
                if x[i] > maxpiv {
                    pivrow = i;
                    maxpiv = x[i];
                }
            }
        }

        // TODO: threshold pivoting.
        let piv = if pivot && !pivots.contains(&pivrow) && maxpiv > d {
            // Swap the max. value in L with the diagonal U(k,k).
            // s.insert(k, RefCell::new(maxpiv));
            x[k] = maxpiv;
            // s.insert(pivrow, RefCell::new(d));

            // Record the pivot in P.
            // row_perm[pivrow] = Some(k);
            // row_perm_inv[k] = Some(pivrow);
            p.swap(pivrow, k);
            pivots.insert(pivrow);

            maxpiv // new diagonal value
        } else {
            d
        };

        // Copy the column elements of U and L, throwing out zeros.
        let u_mat_k: Vec<(usize, f64)> = x
            .iter()
            .enumerate()
            .filter(|(i, _)| *i <= k)
            .filter(|(_, xi)| **xi != 0.0)
            .map(|(i, xi)| (i, *xi))
            .collect();
        let l_mat_k: Vec<(usize, f64)> = x
            .iter()
            .enumerate()
            .filter(|(i, _)| *i > k)
            .filter(|(_, xi)| **xi != 0.0)
            .map(|(i, xi)| (i, *xi))
            .collect();
        u_mat[k] = u_mat_k;
        l_mat[k] = l_mat_k;

        // u_mat[k] = s.range(..k + 1).map(|e| (*e.0, *e.1.borrow())).collect(); // todo: throw out zeros
        // l_mat[k] = s.range(k + 1..).map(|e| (*e.0, *e.1.borrow())).collect();

        // Divide column k of L by U(k,k).
        for l in &mut l_mat[k] {
            l.1 /= piv;
        }
    }

    // Renumber the rows so the data structure represents L and U, not PtL and PtU.
    for row in &mut l_mat {
        for e in row {
            e.0 = p[e.0];
        }
    }

    (l_mat, u_mat, p)
}

fn sp_lsolve_d(l_mat: &Matrix, b: &Col, x: &mut Vec<f64>) {
    // FOR(e, b) x[e->fst] = e->snd;
    // FOR(e, x) FOR(l, L[e->fst])
    //   x[l->fst] -= l->snd * e->snd;

    // > Since we know the nonzero structure of uj before we start,
    // we need only initialize and manipulate the positions in this
    // dense vector that correspond to nonzero positions. TODO
    x.fill(0.0);
    for (bi, bx) in b {
        x[*bi] = *bx; // scatter
    }

    for e0 in 0..x.len() {
        for l in &l_mat[e0] {
            let e1 = x[e0];

            x[l.0] -= l.1 * e1;
        }
    }
}

fn sp_lsolve(l_mat: &Matrix, b: &Col, p: &[usize]) -> BTreeMap<usize, RefCell<f64>> {
    // FOR(e, b) x[e->fst] = e->snd;
    // FOR(e, x) FOR(l, L[e->fst])
    //   x[l->fst] -= l->snd * e->snd;

    let mut x = BTreeMap::<usize, RefCell<f64>>::new();
    for e in b {
        x.insert(e.0, RefCell::new(e.1));
    }

    let mut x_rows = SortedVec::from_unsorted(Vec::from_iter(x.keys().cloned()));

    let mut i = 0;
    while i < x_rows.len() {
        let e0 = x_rows[i];
        // let e0p = p[e0];
        let e0p = e0;

        for l in &l_mat[e0p] {
            let l0p = p[l.0];
            let e1 = x[&e0].borrow().to_owned();

            match x.get(&l0p) {
                Some(xl0_rc) => {
                    let mut xl0 = xl0_rc.try_borrow_mut().unwrap();
                    *xl0 -= l.1 * e1;
                }
                None => {
                    x.insert(l.0, RefCell::new(-l.1 * e1));

                    if l.0 > e0 {
                        x_rows.insert(l.0);
                    }
                }
            };
        }

        i += 1;
    }

    x
}

pub fn lsolve(l_mat: &Matrix, b: &mut [f64]) {
    // FOR(e, b) x[e->fst] = e->snd;
    // FOR(e, x) FOR(l, L[e->fst])
    //   x[l->fst] -= l->snd * e->snd;

    for e0 in 0..b.len() {
        for l in &l_mat[e0] {
            b[l.0] -= l.1 * b[e0];
        }
    }
}

pub fn ltsolve(l_mat: &Matrix, b: &mut [f64]) {
    for e0 in (0..b.len()).rev() {
        for l in l_mat[e0].iter().rev() {
            b[e0] -= l.1 * b[l.0];
        }
    }
}

pub fn usolve(u_mat: &Matrix, b: &mut [f64]) {
    // FORR(e, b) FORR(u, U[e->fst])
    //   if (u->fst == e->fst) e->snd /= u->snd;
    //   else b[u->fst] -= u->snd * e->snd;

    for e0 in (0..b.len()).rev() {
        for u in u_mat[e0].iter().rev() {
            if u.0 == e0 {
                b[e0] /= u.1;
            } else {
                b[u.0] -= u.1 * b[e0];
            }
        }
    }
}

pub fn utsolve(u_mat: &Matrix, b: &mut [f64]) {
    for e0 in 0..b.len() {
        for u in &u_mat[e0] {
            if u.0 == e0 {
                b[e0] /= u.1;
            } else {
                b[e0] -= u.1 * b[u.0];
            }
        }
    }
}
