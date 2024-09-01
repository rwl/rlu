use std::iter::zip;

fn main() {
    let n = 5;
    let mut a_mat: rlu::Matrix<usize, f64> = vec![vec![]; n];
    for i in 0..n {
        if i == n - 1 {
            a_mat[i].push((0, 1.0));
        }
        if i >= 2 {
            a_mat[i].push((i - 2, 3.0));
        }
        if i >= 1 {
            a_mat[i].push((i - 1, 4.0));
        }
        a_mat[i].push((i, 5.0));
        if i + 1 < n {
            a_mat[i].push((i + 1, 2.0));
        }
        if i == 0 {
            a_mat[i].push((n - 1, 1.0));
        }
    }
    // M = |5 4 3     1|
    //     |2 5 4 3    |
    //     |  2 5 4 3  |
    //     |    2 5 4 3|
    //     |      2 5 4|
    //     |1       2 5|
    #[cfg(feature = "debug")]
    print!("A =\n{}", rlu::matrix_table(&a_mat));

    let (n, rowidx, colptr, values) = matrix_to_csc(&a_mat);

    // let p = None;
    let mut q: Vec<usize> = (0..n).map(|i| i).collect();
    // q.swap(0, 0);
    q.swap(0, 1);
    q.swap(2, 3);
    // q.swap(0, 4);
    // q.swap(5, 6);
    println!("q = {:?}", q);

    let b0: Vec<f64> = (1..=n).map(|i| i as f64).collect();
    // x = [1,2,...,n]'
    // let b0: Vec<f64> = (0..n).map(|i| 1.0 + i as f64 / n as f64).collect();
    // x = [1,...,2]'
    println!("b = {:?}", b0);

    // let pivot = false;
    // let pivot = true;

    {
        let mut x = b0.clone();
        //         let (l_mat, u_mat, p) = lu_decomposition(&a_mat, Some(&q), pivot);
        //
        //         println!("p = {:?}", p);
        //
        //         let mut b = vec![0.0; n];
        //         for i in 0..n {
        //             b[p[i].unwrap()] = b0[i];
        //             // b[i] = b0[i];
        //         }
        //
        //         lsolve(&l_mat, &mut b);
        //         usolve(&u_mat, &mut b);
        //
        //         let mut x = vec![0.0; n]; // inverse permutation
        //         for i in 0..n {
        //             x[q[i]] = b[i];
        //             // x[i] = b[i];
        //         }

        rlu::solve(n, &rowidx, &colptr, &values, Some(&q), &mut x, false).unwrap();

        // Matrix-vector multiply b2 = A*x and print residual.
        let mut b2 = vec![0.0; n];
        for j in 0..n {
            for (i, aij) in &a_mat[j] {
                b2[*i] += *aij * x[j];
            }
        }
        b2.iter().for_each(|bi| println!("{}", bi));
        println!(
            "resid: {}",
            zip(b2, &b0)
                .map(|(b2, b0)| f64::abs(b2 - b0))
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap()
        );
    }
    #[cfg(feature = "ignored")]
    {
        let mut b = b0.clone();
        let (l_mat, u_mat): (Matrix, Matrix) = lu_decomposition(&a_mat, None);
        rlu::utsolve(&u_mat, &mut b);
        rlu::ltsolve(&l_mat, &mut b);
        let x = b;

        let mut b2 = vec![0.0; n];
        for i in 0..n {
            for (j, aij) in &a_mat[i] {
                b2[i] += *aij * x[*j];
            }
        }
        println!("{:?}", b2);
    }
}

pub fn matrix_to_csc<I: rlu::Int, S: rlu::Scalar>(
    m: &rlu::Matrix<I, S>,
) -> (usize, Vec<I>, Vec<I>, Vec<S>) {
    let n = m.len();
    let nnz = m.iter().map(|c| c.len()).fold(0, |acc, e| acc + e);
    let mut rowidx: Vec<I> = Vec::with_capacity(nnz);
    let mut colptr: Vec<I> = Vec::with_capacity(n + 1);
    let mut values: Vec<S> = Vec::with_capacity(nnz);

    let mut idxptr: usize = 0;
    for col in m {
        colptr.push(I::from_usize(idxptr));
        for (j, x) in col {
            values.push(*x);
            rowidx.push(*j);
            idxptr += 1
        }
    }
    colptr.push(I::from_usize(idxptr));

    (n, rowidx, colptr, values)
}
