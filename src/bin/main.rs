use lucinda::{lsolve, ltsolve, lu_decomposition, usolve, utsolve, Matrix};

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

    let b: Vec<f64> = (1..=n).map(|i| i as f64).collect();
    // x = [1,2,...,n]'

    {
        let mut b = b.clone();
        let (l_mat, u_mat): (Matrix, Matrix) = lu_decomposition(&a_mat);
        lsolve(&l_mat, &mut b);
        usolve(&u_mat, &mut b);
        let x = b;

        for xi in &x {
            println!("{}", xi);
        }
        let mut b2 = vec![0.0; n];
        for j in 0..n {
            for (i, aij) in &a_mat[j] {
                b2[*i] += *aij * x[j];
            }
        }
        println!("{:?}", b2);
    }
    {
        let mut b = b.clone();
        let (l_mat, u_mat): (Matrix, Matrix) = lu_decomposition(&a_mat);
        utsolve(&u_mat, &mut b);
        ltsolve(&l_mat, &mut b);
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
