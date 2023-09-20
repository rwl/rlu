use crate::{debug, matrix_to_csc, Matrix};

// Depth-first-search workspace.
pub struct DFS {
    root_list: Vec<usize>,
    ptr_list: Vec<usize>,
    flag: Vec<bool>,
}

impl DFS {
    pub fn new(n: usize) -> Self {
        Self {
            root_list: vec![0; n],
            ptr_list: vec![0; n],
            flag: vec![false; n],
        }
    }

    pub fn ludfs(
        &mut self,
        l_mat: &Matrix,
        b_rowidx: &[usize],
        rperm: &Vec<Option<usize>>,
    ) -> &[usize] {
        debug!("b = {:?}", b_rowidx);

        let csgraph = matrix_to_csc(l_mat);
        // debug!("L =\n{}", csgraph.to_table());

        let n = csgraph.cols();
        let indices = csgraph.rowidx();
        let indptr = csgraph.colptr();

        let mut i_rl_start = n;

        self.flag.fill(false);

        for &e0 in b_rowidx {
            // The depth-first search must mark the vertices it
            // has reached, to avoid repeating parts of the search.
            if self.flag[e0] {
                continue;
            }

            self.dfs(e0, indices, indptr, &mut i_rl_start, rperm);
        }
        debug!("found = {:?}", self.root_list[i_rl_start..].to_vec());

        &self.root_list[i_rl_start..]
    }

    // Based on `depth_first_directed` from SciPy v1.11.
    // Modified to support repeated calls with different start nodes
    // and the same permuted input matrix.
    //
    // Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
    // License: BSD, (C) 2012
    pub(crate) fn dfs(
        &mut self,
        head_node: usize,
        indices: &[usize],
        indptr: &[usize],
        i_rl_start: &mut usize,
        rperm: &Vec<Option<usize>>,
    ) {
        // let n = node_list.len();

        // node_list[0] = head_node;
        self.root_list[0] = head_node;
        let mut i_root_opt: Option<usize> = Some(0);
        // let mut i_nl_end = 1;
        // flag[head_node] = true;

        while let Some(mut i_root) = i_root_opt {
            let pnode = self.root_list[i_root];
            let pnode_p = rperm[pnode];

            if !self.flag[pnode] {
                self.flag[pnode] = true;
                match pnode_p {
                    Some(pnode_p) => self.ptr_list[i_root] = indptr[pnode_p],
                    None => self.ptr_list[i_root] = 0,
                }
            }

            let mut no_children = true;

            let indptr1 = self.ptr_list[i_root];
            let indptr2 = match pnode_p {
                Some(pnode_p) => indptr[pnode_p + 1],
                None => 0,
            };
            // debug!(
            //     "pnode = {}, pnode_p = {:?}, p1 = {}, p2 = {}",
            //     pnode, pnode_p, indptr1, indptr2
            // );
            for i in indptr1..indptr2 {
                let cnode = indices[i];
                if self.flag[cnode] {
                    continue;
                } else {
                    self.ptr_list[i_root] = i;
                    i_root += 1;
                    i_root_opt = Some(i_root);
                    self.root_list[i_root] = cnode;
                    // node_list[i_nl_end] = cnode;
                    // flag[cnode] = true;
                    // i_nl_end += 1;

                    // debug!("i_root = {}, cnode = {}", i_root, cnode);
                    no_children = false;
                    break;
                }
            }

            if *i_rl_start == 0 {
                break;
            }

            if no_children {
                i_root_opt = if i_root > 0 { Some(i_root - 1) } else { None };

                *i_rl_start -= 1;
                self.root_list[*i_rl_start] = pnode;
                // debug!("i_rl_start = {}, pnode = {}", *i_rl_start, pnode);
            }
        }
    }
}
