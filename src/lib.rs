//! Sparse LU Decomposition (Gilbert-Peierls)

mod debug;
mod dfs;
mod rlu;
mod traits;

#[cfg(feature = "rayon")]
mod par_solve;

#[cfg(feature = "spsolve")]
pub mod spsolve;

pub use rlu::*;
pub use traits::*;

#[cfg(feature = "rayon")]
pub use par_solve::*;

#[cfg(feature = "debug")]
pub use debug::matrix_table;
