use std::iter::zip;

use lucinda::{lsolve, lu_decomposition, usolve, Matrix};

fn main() {
    let n = 10;
    let mut a_mat: Matrix = vec![vec![]; n];
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

    // let p = None;
    let mut q: Vec<usize> = (0..n).map(|i| i).collect();
    q.swap(0, 1);
    // q.swap(3, 4);
    // q.swap(1, 9);
    // q.swap(5, 6);
    println!("q = {:?}", q);

    let b0: Vec<f64> = (1..=n).map(|i| i as f64).collect();
    // x = [1,2,...,n]'
    // let b0: Vec<f64> = (0..n).map(|i| 1.0 + i as f64 / n as f64).collect();
    // x = [1,...,2]'

    // let pivot = false;
    let pivot = true;

    {
        // let mut b = b0.clone();
        let (l_mat, u_mat, p) = lu_decomposition(&a_mat, Some(&q), pivot);

        println!("p = {:?}", p);

        let mut b = vec![0.0; n];
        for i in 0..n {
            b[p[i]] = b0[i];
            // b[i] = b0[i];
        }

        lsolve(&l_mat, &mut b);
        usolve(&u_mat, &mut b);

        let mut x = vec![0.0; n]; // inverse permutation
        for i in 0..n {
            x[q[i]] = b[i];
        }

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
    {
        let mut b = b0.clone();
        let (l_mat, u_mat): (Matrix, Matrix) = lu_decomposition(&a_mat, None);
        lucinda::utsolve(&u_mat, &mut b);
        lucinda::ltsolve(&l_mat, &mut b);
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
