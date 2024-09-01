use anyhow::{format_err, Result};

use num_traits::NumAssignOps;
use spsolve::Solver;

use crate::{solve, Int, Scalar};

#[derive(Default)]
pub struct RLU {
    pub control: amd::Control,
}

impl<I, S> Solver<I, S> for RLU
where
    I: Int + NumAssignOps,
    S: Scalar,
{
    fn solve(
        &self,
        n: usize,
        a_i: &[I],
        a_p: &[I],
        a_x: &[S],
        b: &mut [S],
        trans: bool,
    ) -> Result<()> {
        let (p, _p_inv, _info) = amd::order::<I>(I::from_usize(n), &a_p, &a_i, &self.control)
            .map_err(|st| format_err!("amd status: {:?}", st))?;

        solve(n, &a_i, &a_p, &a_x, Some(&p), b, trans)
    }
}
