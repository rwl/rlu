use std::collections::BTreeSet;

use crate::{matrix_to_csc, Col, Matrix};

pub fn ludfs(l_mat: &Matrix, b: &Col, rperm: &Vec<Option<usize>>) -> BTreeSet<usize> {
    println!("b = {:?}", b);

    let csgraph = matrix_to_csc(l_mat);
    // print!("L =\n{}", csgraph.to_csr().to_table());
    // let csgraph = csgraph.transpose().to_csr();

    let mut found = BTreeSet::new();
    // let mut found = Vec::new();
    for (bi, _) in b {
        if found.contains(bi) {
            continue;
        }
        let mut node_list = vec![0; csgraph.cols()];
        let length = dfs(
            *bi,
            csgraph.rowidx(),
            csgraph.colptr(),
            &mut node_list,
            rperm,
        );
        for &node in &node_list[..length] {
            // println!("node[{}] = {}", bi, node);
            if !found.contains(&node) {
                // found.insert(0, node);
                found.insert(node);
            }
        }
    }
    // found.sort();
    println!("found = {:?}", found);

    // > The one remaining issue is that the depth-first search must mark the vertices it
    // has reached, to avoid repeating parts of the search.

    found
}

// Based on `depth_first_directed` from SciPy v1.11.
// Modified to support a permuted input matrix.
//
// Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
// License: BSD, (C) 2012
fn dfs(
    head_node: usize,
    indices: &[usize],
    indptr: &[usize],
    node_list: &mut [usize],
    rperm: &Vec<Option<usize>>,
) -> usize {
    let n = node_list.len();

    let mut root_list = vec![0; n];
    let mut flag = vec![false; n];

    node_list[0] = head_node;
    root_list[0] = head_node;
    let mut i_root: isize = 0;
    let mut i_nl_end = 1;
    flag[head_node] = true;

    while i_root >= 0 {
        let pnode = root_list[i_root as usize];
        let pnode_p = rperm[pnode];
        let mut no_children = true;
        let indptr1 = match pnode_p {
            Some(pnode_p) => indptr[pnode_p],
            None => 0,
        };
        let indptr2 = match pnode_p {
            Some(pnode_p) => indptr[pnode_p + 1],
            None => 0,
        };
        for i in indptr1..indptr2 {
            let cnode = indices[i];
            if flag[cnode] {
                continue;
            } else {
                i_root += 1;
                root_list[i_root as usize] = cnode;
                node_list[i_nl_end] = cnode;
                flag[cnode] = true;
                i_nl_end += 1;
                no_children = false;
                break;
            }
        }

        if i_nl_end == n {
            break;
        }

        if no_children {
            i_root -= 1;
        }
    }
    i_nl_end
}
