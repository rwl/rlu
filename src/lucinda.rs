use sparsetools::{csc::CSC, csr::CSR, graph::depth_first_order};

// Simplified Compressed Column Storage
//
// > We shall represent a column vector as a sequence of records, each containing a
// value and a row index. The row indices need not be in increasing order. We shall
// represent a matrix as an array of column vectors indexed from 0 to n.
pub type Record = (usize, f64);
pub type Col = Vec<Record>;
pub type Matrix = Vec<Col>;

pub fn matrix_to_csr(m: &Matrix) -> CSR<usize, f64> {
    let n = m.len();
    let mut rowidx = vec![];
    let mut colptr = vec![]; // n + 1
    let mut data = vec![];

    let mut idxptr: usize = 0;
    for col in m {
        colptr.push(idxptr);
        for (j, x) in col {
            data.push(*x);
            rowidx.push(*j);
            idxptr += 1
        }
    }
    colptr.push(idxptr);

    let csc = CSC::new(n, n, rowidx, colptr, data).unwrap();

    csc.to_csr()
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
pub fn lu_decomposition(
    a_mat: &Matrix,
    col_perm: Option<&[usize]>,
    pivot: bool,
) -> (Matrix, Matrix, Vec<Option<usize>>) {
    let n = a_mat.len();

    // row_perm(r) = Some(s) means row r of A is row s < jcol of PA (LU = PA).
    // row_perm(r) = None means row r of A has not yet been used as a
    // pivot and is therefore still below the diagonal.
    let mut row_perm: Vec<Option<usize>> = vec![None; n];

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
        println!("\nk = {}, kp = {}", k, kp);

        // Depth-first search from each above-diagonal nonzero of column jcol of A.
        let found = ludfs(&l_mat, &a_mat[kp]);

        // Compute the values of column jcol of L and U in the *dense* vector,
        // allocating storage for fill in L as necessary.
        lucomp(&l_mat, &a_mat[kp], &mut x, &row_perm, &found); // s = L\A[:,j]

        let d = x[k]; // TODO: numerically zero diagonal element at column check

        // Partial pivoting, diagonal elt. has max. magnitude in L.
        let mut pivt = found
            .iter()
            .filter(|i| row_perm[**i].is_none())
            .map(|&i| (i, x[i].abs()))
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .expect("pivot must exist");

        println!("pivrow = {}, maxpiv = {:.6}, d = {:.6}", pivt.0, pivt.1, d);

        // TODO: Threshold pivoting.
        if !pivot || (row_perm[k].is_none() && d.abs() >= pivt.1) {
            // No pivoting, diagonal element has irow = jcol.
            pivt = (k, d);
        };

        // Copy the column elements of U, throwing out zeros.
        for &i in &found {
            if let Some(ip) = row_perm[i] {
                u_mat[k].push((ip, x[i]));
            }
        }

        // Swap the max. value in L with the diagonal U(k,k).
        u_mat[k].push((k, pivt.1));

        // Record the pivot in P.
        row_perm[pivt.0] = Some(k);

        // Copy the column elements of L, throwing out zeros.
        for &i in &found {
            if row_perm[i].is_none() {
                // Divide column k of L by U(k,k).
                l_mat[k].push((i, x[i] / pivt.1));
            }
        }

        // println!("U[:,k] = {:?}", u_mat[k]);
        // println!("L[:,k] = {:?}", l_mat[k]);

        // {
        //     let csgraph = matrix_to_csr(&u_mat);
        //     print!("U =\n{}", csgraph.to_table());
        // }
        // {
        //     let csgraph = matrix_to_csr(&l_mat);
        //     print!("L =\n{}", csgraph.to_table());
        // }
    }

    // Renumber the rows so the data structure represents L and U, not PtL and PtU.
    for row in &mut l_mat {
        for e in row {
            match row_perm[e.0] {
                Some(e0) => e.0 = e0,
                None => panic!(),
            }
        }
    }
    println!("L =\n{}", matrix_to_csr(&l_mat).to_table());

    (l_mat, u_mat, row_perm)
}

fn ludfs(l_mat: &Matrix, b: &Col) -> Vec<usize> {
    println!("b = {:?}", b);

    let csgraph = matrix_to_csr(l_mat);
    print!("L =\n{}", csgraph.to_table());

    // let mut found = HashSet::new();
    let mut found = Vec::new();
    for (bi, _) in b {
        let (nodes, _) = depth_first_order(&csgraph, *bi, false).unwrap();
        for node in nodes {
            // println!("node[{}] = {}", bi, node);
            if !found.contains(&node) {
                // found.insert(0, node);
                found.push(node);
            }
        }
    }
    found.sort();
    println!("found = {:?}", found);

    // > The one remaining issue is that the depth-first search must mark the vertices it
    // has reached, to avoid repeating parts of the search.

    found
}

fn lucomp(
    l_mat: &Matrix,
    b: &Col,
    x: &mut Vec<f64>,
    rperm: &Vec<Option<usize>>,
    found: &Vec<usize>,
) {
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

    // for e0 in 0..x.len() {
    for j in found {
        let e0 = match rperm[*j] {
            Some(jp) => jp,
            None => continue,
        };
        for l in &l_mat[e0] {
            let e1 = x[e0];

            x[l.0] -= l.1 * e1;
        }
    }
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
