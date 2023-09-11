use std::{cell::RefCell, collections::BTreeMap};

use sorted_vec::SortedVec;

// Simplified Compressed Column Storage
//
// > We shall represent a column vector as a sequence of records, each containing a
// value and a row index. The row indices need not be in increasing order. We shall
// represent a matrix as an array of column vectors indexed from 0 to n.
pub type Record = (usize, f64);
pub type Col = Vec<Record>;
pub type Matrix = Vec<Col>;

// LU decomposition without reordering (Gilbert-Peierls)
//
// Note: A is LU-decomposable <=> all principal minors are nonsingular
fn sp_lsolve(l_mat: &Matrix, b: &Col) -> BTreeMap<usize, RefCell<f64>> {
    // FOR(e, b) x[e->fst] = e->snd;
    // FOR(e, x) FOR(l, L[e->fst])
    //   x[l->fst] -= l->snd * e->snd;

    let mut x = BTreeMap::<usize, RefCell<f64>>::new();
    for e in b {
        x.insert(e.0, RefCell::new(e.1));
    }

    let mut x_cols = SortedVec::from_unsorted(Vec::from_iter(x.keys().cloned()));

    let mut i = 0;
    while i < x_cols.len() {
        let e0 = x_cols[i];

        for l in &l_mat[e0] {
            let e1 = x[&e0].borrow().to_owned();

            match x.get(&l.0) {
                Some(xl0_rc) => {
                    let mut xl0 = xl0_rc.try_borrow_mut().unwrap();
                    *xl0 -= l.1 * e1;
                }
                None => {
                    x.insert(l.0, RefCell::new(-l.1 * e1));

                    if l.0 > e0 {
                        x_cols.insert(l.0);
                    }
                }
            };
        }

        i += 1;
    }

    x
}

pub fn lsolve(l_mat: &Matrix, b: &mut Vec<f64>) {
    // FOR(e, b) x[e->fst] = e->snd;
    // FOR(e, x) FOR(l, L[e->fst])
    //   x[l->fst] -= l->snd * e->snd;

    for e0 in 0..b.len() {
        for l in &l_mat[e0] {
            b[l.0] -= l.1 * b[e0];
        }
    }
}

pub fn usolve(u_mat: &Matrix, b: &mut Vec<f64>) {
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

// 1. for j:= to n do
// 2.   {Compute column j of U and L.}
// 3.   Solve Ljuj = aj for uj;
// 4.   b'j := a'j - L'juj;
// 5.   Pivot: Swap bjj with the largest-magnitude element of b'j;
// 6.   ujj := bjj;
// 7.   l'j := b'j / ujj;
// 8. od

pub fn lu_decomposition(a_mat: &Matrix) -> (Matrix, Matrix) {
    let n = a_mat.len();
    let mut l_mat: Matrix = vec![vec![]; n];
    let mut u_mat: Matrix = vec![vec![]; n];
    for j in 0..n {
        let s: BTreeMap<usize, RefCell<f64>> = sp_lsolve(&l_mat, &a_mat[j]);

        // No pivoting, diagonal element has irow = jcol.
        // Partial pivoting, diagonal elt. has max. magnitude in L.

        let d = match s.get(&j) {
            Some(d) => d.borrow(),
            None => s.last_key_value().unwrap().1.borrow(), // numerically zero diagonal element at column
        };

        // Copy the column elements of U and L, throwing out zeros.
        u_mat[j] = s.range(..j + 1).map(|e| (*e.0, *e.1.borrow())).collect(); // todo: throw out zeros
        l_mat[j] = s.range(j + 1..).map(|e| (*e.0, *e.1.borrow())).collect();

        // todo: Diagonal element has been found. Swap U(jcol,jcol) from L into U.
        // Record the pivot in P.

        // Divide column j of L by U(j,j).
        for l in &mut l_mat[j] {
            l.1 /= *d;
        }
    }
    (l_mat, u_mat)
}
