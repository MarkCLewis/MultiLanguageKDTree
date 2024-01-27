use std::{fs::File, io::Write};

use core_simd::simd::*;

use crate::simd_particle::*;

use super::simd_particle::Particle;

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
        cm: f64x4,
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

pub fn allocate_node_vec(num_parts: usize) -> Vec<KDTree> {
    let num_nodes = 2 * (num_parts / (MAX_PARTS - 1) + 1);
    let mut ret = Vec::new();
    ret.resize(num_nodes, KDTree::leaf(0, NEGS));
    ret
}

// Returns the index of the last Node used in the construction.
pub fn build_tree(
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
        let mut min = f64x4::splat(1e100);
        let mut max = f64x4::splat(-1e100);
        let mut m = 0.0;
        let mut cm = f64x4::splat(0.0);
        for i in start..end {
            m += particles[indices[i]].m;
            cm += f64x4::splat(particles[indices[i]].m) * particles[indices[i]].p;
            min = f64x4::simd_min(min, particles[indices[i]].p);
            max = f64x4::simd_max(max, particles[indices[i]].p);
        }
        cm /= f64x4::splat(m);
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
            use std::cmp::Ordering::*;
            match high.cmp(&mid) {
                Less => s = high + 1,
                Greater => e = high,
                Equal => s = e,
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

fn accel_recur(cur_node: usize, p: usize, particles: &Vec<Particle>, nodes: &Vec<KDTree>) -> f64x4 {
    // println!("accel {}", cur_node);
    match nodes[cur_node] {
        KDTree::Leaf { num_parts, leaf_parts} => {
            let mut acc = f64x4::splat(0.0);
            for i in 0..num_parts {
                if leaf_parts[i] != p {
                    acc += calc_pp_accel(&particles[p], &particles[leaf_parts[i]]);
                }
            }
            acc
        }
        KDTree::Internal { m, cm, size, left, right, .. } => {
            let dp = particles[p].p - cm;
            let dp2 = dp * dp;
            let dist_sqr = dp2.reduce_sum();
            // println!("dist = {}, size = {}", dist, size);
            if size * size < THETA * THETA * dist_sqr {
                let dist = f64::sqrt(dist_sqr);
                let magi = f64x4::splat(-m / (dist_sqr * dist));
                dp * magi
            } else {
                accel_recur(left, p, particles, nodes)
                    + accel_recur(right, p, particles, nodes)
            }
        }
    }
}

pub fn calc_accel(p: usize, particles: &Vec<Particle>, nodes: &Vec<KDTree>) -> f64x4 {
    accel_recur(0, p, particles, nodes)
}

pub fn simple_sim(bodies: &mut Vec<Particle>, dt: f64, steps: i64) {
    let dt_vec = f64x4::splat(dt);
    let mut acc = Vec::new();
    for _ in 0..bodies.len() {
        acc.push(f64x4::splat(0.0))
    }
    // let mut time = Instant::now();
    let mut tree = allocate_node_vec(bodies.len());
    let mut indices: Vec<usize> = (0..bodies.len()).collect();

    for _ in 0..steps {
        // if step % 100 == 0 {
        //     let elapsed_secs = time.elapsed().as_nanos() as f64 / 1e9;
        //     println!("Step = {}, duration = {}, n = {}, nodes = {}", step, elapsed_secs, bodies.len(), tree.len());
        //     time = Instant::now();
        // }
        for (new_ix, old_ix) in indices.iter_mut().enumerate().take(bodies.len()) {
            *old_ix = new_ix;
        }

        build_tree(&mut indices, 0, bodies.len(), bodies, 0, &mut tree);
        // if step % 10 == 0 {
        //     print_tree(step, &tree, &bodies);
        // }
        for (index, item) in acc.iter_mut().enumerate().take(bodies.len()) {
            *item = calc_accel(index, bodies, &tree)
        }

        for i in 0..bodies.len() {
            bodies[i].v += dt_vec * acc[i];
            let dp = dt_vec * bodies[i].v;
            bodies[i].p += dp;
            acc[i] = f64x4::splat(0.0);
        }
    }
}

#[allow(dead_code)]
fn print_tree(step: i64, tree: &Vec<KDTree>, particles: &[Particle]) -> std::io::Result<()> {
    let mut file = File::create(format!("tree{step}.txt"))?;

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

#[allow(dead_code)]
fn recur_test_tree_struct(
    node: usize,
    nodes: &Vec<KDTree>,
    particles: &Vec<Particle>,
    mut min: f64x4,
    mut max: f64x4,
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
    use crate::{kd_tree, simd_particle};
    use core_simd::simd::*;

    #[test]
    fn single_node() {
        let parts = simd_particle::two_bodies();
        let mut node_vec = kd_tree::allocate_node_vec(parts.len());
        assert_eq!(node_vec.len(), 2);
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
        match node_vec[0] {
            kd_tree::KDTree::Leaf { num_parts, .. } => assert_eq!(num_parts, parts.len()),
            _ => assert!(false, "Root isn't leaf of right size when small.")
        };
    }

    #[test]
    fn two_leaves() {
        let parts = simd_particle::circular_orbits(11);
        let mut node_vec = kd_tree::allocate_node_vec(parts.len());
        assert_eq!(node_vec.len(), 6);
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
        kd_tree::recur_test_tree_struct(
            0,
            &node_vec,
            &parts,
            f64x4::splat(-1e100),
            f64x4::splat(1e100),
        );
        assert!(std::matches!(node_vec[0], kd_tree::KDTree::Internal { .. }));
        match (node_vec[1], node_vec[2]) {
            (kd_tree::KDTree::Leaf { num_parts: n1, ..}, kd_tree::KDTree::Leaf {num_parts: n2, ..}) => {
                assert_eq!(n1 + n2, 12);
            }
            _ => assert!(false, "Node vectors weren't leaves.")
        }
    }

    #[test]
    fn big_solar() {
        let parts = simd_particle::circular_orbits(5000);
        let mut node_vec = kd_tree::allocate_node_vec(parts.len());
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
        kd_tree::recur_test_tree_struct(
            0,
            &node_vec,
            &parts,
            f64x4::splat(-1e100),
            f64x4::splat(1e100),
        );
    }

    #[test]
    fn big_solar_with_steps() {
        let mut parts = simd_particle::circular_orbits(5000);
        kd_tree::simple_sim(&mut parts, 1e-3, 10);

        let mut node_vec = kd_tree::allocate_node_vec(parts.len());
        let mut indices: Vec<usize> = (0..parts.len()).collect();
        kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
        kd_tree::recur_test_tree_struct(
            0,
            &node_vec,
            &parts,
            f64x4::splat(-1e100),
            f64x4::splat(1e100),
        );
    }
}
