fn main() {
    let n = 4;

    // int n = 4;
    // int Ap[5] = {0, 1, 3, 6, 8}; // Column pointers
    // int Ai[8] = {0, 1, 2, 1, 2, 3, 2, 3}; // Row indices
    // double Ax[8] = {10, 2, 20, 3, 30, 4, 5, 40}; // Non-zero values

    let a_colptr = vec![0, 1, 3, 6, 8];
    let a_rowidx = vec![0, 1, 2, 1, 2, 3, 2, 3];
    let a_values = vec![10.0, 2.0, 20.0, 3.0, 30.0, 4.0, 5.0, 40.0];
    let b = vec![10.0, 40.0, 90.0, 160.0];

    let (p, searches, rp_inv) = {
        // let mut x = b.clone();
        // rlu::solve::<usize, f64, usize>(n, &a_rowidx, &a_colptr, &a_values, None, &mut x, false);

        let mut founds = Vec::with_capacity(n);
        let mut rp_inv = vec![0; n];

        let (l_mat, u_mat, rp) = rlu::lu_decomposition::<usize, f64, usize>(
            n,
            &a_rowidx,
            &a_colptr,
            &a_values,
            None,
            Some(&mut founds),
            Some(&mut rp_inv),
            true,
        );

        let mut x = vec![0.0; n];
        for i in 0..n {
            x[rp[i].unwrap()] = b[i];
        }

        rlu::lsolve(&l_mat, &mut x);
        rlu::usolve(&u_mat, &mut x);

        // b.copy_from_slice(&x);

        // [1.0, -970.0, 660.0, -62.0]
        println!("X = {:?}", x);

        (rp, founds, rp_inv)
    };

    println!("FOUND = {:?}", searches);

    {
        // let a_values2 = vec![15.0, 1.0, 25.0, 4.0, 35.0, 6.0, 2.0, 45.0];

        let (l_mat, u_mat) = rlu::lu_redecomposition::<usize, f64, usize>(
            n,
            &a_rowidx,
            &a_colptr,
            &a_values,
            &searches, /*, &p*/
            Some(&rp_inv),
            None,
            true,
        );

        let mut x = vec![0.0; n];
        for i in 0..n {
            x[p[i].unwrap()] = b[i];
        }

        rlu::lsolve(&l_mat, &mut x);
        rlu::usolve(&u_mat, &mut x);

        // b.copy_from_slice(&x);

        println!("X2 = {:?}", x);
    }
}
