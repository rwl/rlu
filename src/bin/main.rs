use lucinda::{lsolve, lu_decomposition, usolve, Matrix, Row};

fn main() {
    let n = 100;
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

    let mut x: Row = vec![];
    for i in 0..n {
        x.push((i, i as f64 + 1.0));
    }
    // x = [1,2,...,n]'

    let (l_mat, u_mat): (Matrix, Matrix) = lu_decomposition(&a_mat);
    x = usolve(&u_mat, &mut lsolve(&l_mat, &x));

    for e in x {
        println!("{} {}", e.0, e.1);
    }
}
