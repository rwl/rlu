use crate::traits::{Int, Scalar};
use crate::{lsolve, ltsolve, lu_decomposition, usolve, utsolve};

use anyhow::{format_err, Result};
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSliceMut;

/// Solve `Ax=b` for one or more right-hand-sides given the numeric
/// factorization of A from `factor`.
pub fn par_solve<I: Int + Sync, S: Scalar + Send + Sync, P: Int + Sync>(
    n: usize,
    a_rowidx: &[I],
    a_colptr: &[I],
    a_values: &[S],
    col_perm: Option<&[P]>,
    rhs: &mut [S],
    trans: bool,
) -> Result<()> {
    if rhs.len() % n != 0 {
        return Err(format_err!(
            "len rhs ({}) must be a multiple of n ({})",
            rhs.len(),
            n
        ));
    }

    let (l_mat, u_mat, p) =
        lu_decomposition(n, a_rowidx, a_colptr, a_values, col_perm, None, None, true);

    rhs.par_chunks_exact_mut(n)
        .try_for_each_with(vec![S::zero(); n], |x, b| -> Result<()> {
            for i in 0..n {
                x[p[i].unwrap()] = b[i];
            }

            if !trans {
                lsolve(&l_mat, x);
                usolve(&u_mat, x);
            } else {
                utsolve(&u_mat, x);
                ltsolve(&l_mat, x);
            }

            match col_perm {
                Some(cperm) => {
                    for i in 0..n {
                        b[cperm[i].to_index()] = x[i];
                    }
                }
                None => b.copy_from_slice(&x),
            }

            Ok(())
        })
}
