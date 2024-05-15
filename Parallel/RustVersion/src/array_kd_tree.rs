use std::{fs::File, io::Write};

use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
use rayon::prelude::*;

use crate::array_particle::*;

//use super::array_particle::Particle;

const MAX_PARTS: usize = 7;
const THETA: f64 = 0.3;
const NEGS: [usize; MAX_PARTS] = [usize::MAX; MAX_PARTS];

#[derive(Clone, Copy)]
pub enum KDTree {
    Leaf {
        num_parts: usize,
        leaf_parts: [usize; MAX_PARTS]
    },

    Internal {
        split_dim: usize,
        split_val: f64,
        m: f64,
        cm: [f64; 3],
        size: f64,
        left: usize,
        right: usize
    }
}

impl KDTree {
    pub fn leaf<'a>(num_parts: usize, particles: [usize; MAX_PARTS]) -> KDTree {
        KDTree::Leaf {
            num_parts: num_parts,
            leaf_parts: particles,
        }
    }
}

fn nodes_needed_for_particles(num_parts: usize) -> usize {
    2 * (num_parts / (MAX_PARTS / 2) + 1)
}

pub fn allocate_node_vec(num_parts: usize) -> Vec<KDTree> {
    let num_nodes = nodes_needed_for_particles(num_parts);
    let mut ret = Vec::new();
    ret.resize(num_nodes, KDTree::leaf(0, NEGS));
    ret
}

// Returns the index of the last Node used in the construction.
pub fn build_tree<'a>(
    indices: &mut Vec<usize>,
    start: usize,
    end: usize,
    particles: &Vec<Particle>,
    cur_node: usize,
    nodes: &mut Vec<KDTree>,
) -> usize {
    // println!("start = {} end = {} cur_node = {}", start, end, cur_node);
    let np = end - start;
    // println!("s = {}, e = {}, cn = {}", start, end, cur_node);
    if np <= MAX_PARTS {
        if cur_node >= nodes.len() {
            nodes.resize(cur_node + 1, KDTree::leaf(0, NEGS));
        }
        let mut parts = [0; MAX_PARTS];
        for i in 0..np {
            parts[i] = indices[start + i]
        }
        nodes[cur_node] = KDTree::Leaf { num_parts: np, leaf_parts: parts };
        cur_node
    } else {
        // Pick split dim and value
        let mut min = [1e100, 1e100, 1e100];
        let mut max = [-1e100, -1e100, -1e100];
        let mut m = 0.0;
        let mut cm = [0.0, 0.0, 0.0];
        for i in start..end {
            m += particles[indices[i]].m;
            cm[0] += particles[indices[i]].m * particles[indices[i]].p[0];
            cm[1] += particles[indices[i]].m * particles[indices[i]].p[1];
            cm[2] += particles[indices[i]].m * particles[indices[i]].p[2];
            min[0] = f64::min(min[0], particles[indices[i]].p[0]);
            min[1] = f64::min(min[1], particles[indices[i]].p[1]);
            min[2] = f64::min(min[2], particles[indices[i]].p[2]);
            max[0] = f64::max(max[0], particles[indices[i]].p[0]);
            max[1] = f64::max(max[1], particles[indices[i]].p[1]);
            max[2] = f64::max(max[2], particles[indices[i]].p[2]);
        }
        cm[0] /= m;
        cm[1] /= m;
        cm[2] /= m;
        let mut split_dim = 0;
        for dim in 1..3 {
            if max[dim] - min[dim] > max[split_dim] - min[split_dim] {
                split_dim = dim
            }
        }
        let size = max[split_dim] - min[split_dim];

        // Partition particles on split_dim
        let mid = (start + end) / 2;
        let mut s = start;
        let mut e = end;
        while s + 1 < e {
            let pivot = fastrand::usize(s..e);
            indices.swap(s, pivot);
            let mut low = s + 1;
            let mut high = e - 1;
            while low <= high {
                if particles[indices[low]].p[split_dim] < particles[indices[s]].p[split_dim] {
                    low += 1;
                } else {
                    indices.swap(low, high);
                    high -= 1;
                }
            }
            indices.swap(s, high);
            if high < mid {
                s = high + 1;
            } else if high > mid {
                e = high;
            } else {
                s = e;
            }
        }
        let split_val = particles[indices[mid]].p[split_dim];

        // Recurse on children and build this node.
        let left = build_tree(indices, start, mid, particles, cur_node + 1, nodes);
        let right = build_tree(indices, mid, end, particles, left + 1, nodes);

        if cur_node >= nodes.len() {
            nodes.resize(cur_node + 1, KDTree::leaf(0, NEGS));
        }
        nodes[cur_node] = KDTree::Internal { split_dim, split_val, m, cm, size, left: cur_node + 1, right: left+1 };

        right
    }
}

fn accel_recur(cur_node: usize, p: usize, particles: &Vec<Particle>, nodes: &Vec<KDTree>) -> [f64; 3] {
    // println!("accel {}", cur_node);
    match nodes[cur_node] {
        KDTree::Leaf { num_parts, leaf_parts} => {
            let mut acc = [0.0, 0.0, 0.0];
            for i in 0..(num_parts) {
                if leaf_parts[i] != p {
                    let pp_acc = calc_pp_accel(&particles[p], &particles[leaf_parts[i]]);
                    acc[0] += pp_acc[0];
                    acc[1] += pp_acc[1];
                    acc[2] += pp_acc[2];
                }
            }
            acc
        }
        KDTree::Internal { m, cm, size, left, right, .. } => {
            let dx = particles[p].p[0] - cm[0];
            let dy = particles[p].p[1] - cm[1];
            let dz = particles[p].p[2] - cm[2];
            let dist_sqr = dx * dx + dy * dy + dz * dz;
            // println!("dist = {}, size = {}", dist, nodes[cur_node].size);
            if size * size < THETA * THETA * dist_sqr {
                let dist = f64::sqrt(dist_sqr);
                let magi = -m / (dist_sqr * dist);
                [dx * magi, dy * magi, dz * magi]
            } else {
                let left_acc = accel_recur(left, p, particles, nodes);
                let right_acc = accel_recur(right, p, particles, nodes);
                [left_acc[0] + right_acc[0], left_acc[1] + right_acc[1], left_acc[2] + right_acc[2]]
            }
        }
    }
}

pub fn calc_accel(p: usize, particles: &Vec<Particle>, nodes: &Vec<KDTree>) -> [f64; 3] {
    accel_recur(0, p, particles, nodes)
}

pub fn simple_sim(bodies: &mut Vec<Particle>, dt: f64, steps: i64) {
    let mut acc = Vec::new();
    for _ in 0..bodies.len() {
        acc.push([0.0, 0.0, 0.0])
    }
    // let mut time = Instant::now();
    let mut tree = allocate_node_vec(bodies.len());
    let mut indices: Vec<usize> = (0..bodies.len()).collect();

    for _step in 0..steps {
        // if step % 100 == 0 {
        //     let elapsed_secs = time.elapsed().as_nanos() as f64 / 1e9;
        //     println!("Step = {}, duration = {}, n = {}, nodes = {}", step, elapsed_secs, bodies.len(), tree.len());
        //     time = Instant::now();
        // }
        // for i in 0..bodies.len() {
        //     indices[i] = i;
        // }
        indices.par_iter_mut().enumerate().for_each(|(i, ind)| *ind = i);
        build_tree(&mut indices, 0, bodies.len(), bodies, 0, &mut tree);
        // if step % 10 == 0 {
        //     print_tree(step, &tree, &bodies);
        // }
        acc.par_iter_mut().enumerate().for_each(|(i, acc)| *acc = calc_accel(i, &bodies, &tree));

        bodies.par_iter_mut().zip(&mut acc).for_each(|(b, a)| {
            b.v[0] += dt * a[0];
            b.v[1] += dt * a[1];
            b.v[2] += dt * a[2];
            let dx = dt * b.v[0];
            let dy = dt * b.v[1];
            let dz = dt * b.v[2];
            b.p[0] += dx;
            b.p[1] += dy;
            b.p[2] += dz;
            a[0] = 0.0;
            a[1] = 0.0;
            a[2] = 0.0;
        });
    }
}

fn print_tree(step: i64, tree: &Vec<KDTree>, particles: &Vec<Particle>) -> std::io::Result<()> {
    let mut file = File::create(format!("tree{}.txt", step))?;

    file.write_fmt(format_args!("{}\n", tree.len()))?;
    for n in tree {
        match n {
            KDTree::Leaf { num_parts, leaf_parts } => {
                file.write_fmt(format_args!("L {}\n", num_parts))?;
                for i in 0..*num_parts {
                    let p = leaf_parts[i];
                    file.write_fmt(format_args!(
                        "{} {} {}\n",
                        particles[p].p[0], particles[p].p[1], particles[p].p[2]
                    ))?;
                }
            }
            KDTree::Internal { split_dim, split_val, left, right, .. } => {
                file.write_fmt(format_args!(
                    "I {} {} {} {}\n",
                    split_dim, split_val, left, right
                ))?;
            }
        }
    }

    Ok(())
}

fn recur_test_tree_struct(
    node: usize,
    nodes: &Vec<KDTree>,
    particles: &Vec<Particle>,
    mut min: [f64; 3],
    mut max: [f64; 3],
) {
    match nodes[node] {
        KDTree::Leaf { num_parts, leaf_parts } => {
            for index in 0..num_parts {
                let i = leaf_parts[index];
                for dim in 0..2 {
                    assert!(
                        particles[i].p[dim] >= min[dim],
                        "Particle dim {} is below min. i={} p={} min={}",
                        dim,
                        i,
                        particles[i].p[dim],
                        min[dim]
                    );
                    assert!(
                        particles[i].p[dim] < max[dim],
                        "Particle dim {} is above max. i={} p={} max={}",
                        dim,
                        i,
                        particles[i].p[dim],
                        max[dim]
                    );
                }
            }
        }
        KDTree::Internal { split_dim, split_val, left, right, .. } => {
            let split_dim = split_dim;
            let tmin = min[split_dim];
            let tmax = max[split_dim];
            max[split_dim] = split_val;
            recur_test_tree_struct(left, nodes, particles, min, max);
            max[split_dim] = tmax;
            min[split_dim] = split_val;
            recur_test_tree_struct(right, nodes, particles, min, max);
            min[split_dim] = tmin;
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{array_kd_tree, array_particle};

    #[test]
    fn single_node() {
        let parts = array_particle::two_bodies();
        let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
        assert_eq!(node_vec.len(), 2);
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        array_kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
        match node_vec[0] {
            array_kd_tree::KDTree::Leaf { num_parts, .. } => assert_eq!(num_parts, parts.len()),
            _ => assert!(false, "Root isn't leaf of right size when small.")
        };
    }

    #[test]
    fn two_leaves() {
        let parts = array_particle::circular_orbits(11);
        let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
        assert_eq!(node_vec.len(), 6);
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        array_kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
        array_kd_tree::recur_test_tree_struct(
            0,
            &node_vec,
            &parts,
            [-1e100, -1e100, -1e100],
            [1e100, 1e100, 1e100],
        );
        assert!(std::matches!(node_vec[0], array_kd_tree::KDTree::Internal { .. }));
        match (node_vec[1], node_vec[2]) {
            (array_kd_tree::KDTree::Leaf { num_parts: n1, ..}, array_kd_tree::KDTree::Leaf {num_parts: n2, ..}) => {
                assert_eq!(n1 + n2, 12);
            }
            _ => assert!(false, "Node vectors weren't leaves.")
        }
    }

    #[test]
    fn big_solar() {
        let parts = array_particle::circular_orbits(5000);
        let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        array_kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
        array_kd_tree::recur_test_tree_struct(
            0,
            &node_vec,
            &parts,
            [-1e100, -1e100, -1e100],
            [1e100, 1e100, 1e100],
        );
    }

    #[test]
    fn big_solar_with_steps() {
        let mut parts = array_particle::circular_orbits(5000);
        array_kd_tree::simple_sim(&mut parts, 1e-3, 10);

        let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        array_kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
        array_kd_tree::recur_test_tree_struct(
            0,
            &node_vec,
            &parts,
            [-1e100, -1e100, -1e100],
            [1e100, 1e100, 1e100],
        );
    }
}
