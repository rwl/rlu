use std::{cell::RefCell, cmp::Reverse, collections::BTreeMap};

use sorted_vec::{ReverseSortedVec, SortedVec};

// Simplified Compressed Row Storage
pub type Entry = (usize, f64);
pub type Row = Vec<Entry>;
pub type Matrix = Vec<Row>;

// LU decomposition without reordering (Gilbert-Peierls)
//
// Note: A is LU-decomposable <=> all principal minors are nonsingular
pub fn lsolve(l_mat: &Matrix, b: &Row) -> BTreeMap<usize, RefCell<f64>> {
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
                    } else {
                        println!("{} {}", l.0, e0);
                    }
                }
            };
        }

        i += 1;
    }

    x
}

pub fn usolve(u_mat: &Matrix, b: &mut BTreeMap<usize, RefCell<f64>>) -> Row {
    // FORR(e, b) FORR(u, U[e->fst])
    //   if (u->fst == e->fst) e->snd /= u->snd;
    //   else b[u->fst] -= u->snd * e->snd;

    let mut b_cols =
        ReverseSortedVec::from_unsorted(Vec::from_iter(b.keys().into_iter().map(|&k| Reverse(k))));

    let mut i = 0;
    while i < b_cols.len() {
        let e0 = b_cols[i];

        for u in u_mat[e0.0].iter().rev() {
            if u.0 == e0.0 {
                let e1_rc = b.get_mut(&e0.0).unwrap();
                let mut e1 = e1_rc.borrow_mut();
                *e1 /= u.1;
            } else {
                let e1 = b[&e0.0].borrow().to_owned();

                match b.get_mut(&u.0) {
                    Some(bu0_rc) => {
                        let mut bu0 = bu0_rc.borrow_mut();
                        *bu0 -= u.1 * e1;
                    }
                    None => {
                        b.insert(u.0, RefCell::new(-u.1 * e1));

                        if u.0 > e0.0 {
                            b_cols.insert(Reverse(u.0));
                        } else {
                            println!("{} {}", u.0, e0.0);
                        }
                    }
                }
            }
        }

        i += 1;
    }
    let mut x: Row = Vec::with_capacity(b.len());
    for e in b {
        x.push((*e.0, e.1.borrow().to_owned()));
    }
    x
}

pub fn lu_decomposition(a_mat: &Matrix) -> (Matrix, Matrix) {
    let n = a_mat.len();
    let mut l_mat: Matrix = vec![vec![]; n];
    let mut u_mat: Matrix = vec![vec![]; n];
    for k in 0..n {
        let s: BTreeMap<usize, RefCell<f64>> = lsolve(&l_mat, &a_mat[k]);
        let d = match s.get(&k) {
            Some(d) => d.borrow(),
            None => s.last_key_value().unwrap().1.borrow(),
        };
        u_mat[k] = s.range(..k + 1).map(|e| (*e.0, *e.1.borrow())).collect();
        l_mat[k] = s.range(k + 1..).map(|e| (*e.0, *e.1.borrow())).collect();

        for l in &mut l_mat[k] {
            l.1 /= *d;
        }
    }
    (l_mat, u_mat)
}
