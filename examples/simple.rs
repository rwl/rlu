fn main() {
    let n = 3;

    let colptr = vec![0, 2, 5, 6];
    let rowind = vec![0, 1, 0, 1, 2, 2];
    let data = vec![4.0, 1.0, 1.0, 3.0, 1.0, 2.0];
    let mut b = vec![6.0, 7.0, 8.0];
    rlu::solve::<usize, f64, usize>(n, &rowind, &colptr, &data, None, &mut b, false);

    println!("{:?}", b);

    let rowptr = vec![0, 2, 4, 6];
    let colind = vec![0, 1, 0, 1, 1, 2];
    let data = vec![4.0, 1.0, 1.0, 3.0, 1.0, 2.0];
    let mut b = vec![6.0, 7.0, 8.0];
    rlu::solve::<usize, f64, usize>(n, &colind, &rowptr, &data, None, &mut b, true);

    println!("{:?}", b);
}
