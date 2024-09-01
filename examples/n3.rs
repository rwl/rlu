fn main() {
    let n = 3;

    let rowptr = vec![0, 2, 4, 6];
    let colind = vec![0, 1, 0, 1, 1, 2];
    let data = vec![4.0, 1.0, 1.0, 3.0, 1.0, 2.0];

    let b = vec![6.0, 7.0, 8.0];
    let mut rhs = [b.clone(), b.clone()].concat();

    rlu::par_solve::<usize, f64, usize>(n, &colind, &rowptr, &data, None, &mut rhs, true).unwrap();

    rhs.chunks_exact(n).for_each(|x| {
        println!("{:?}", x);
    })
}
