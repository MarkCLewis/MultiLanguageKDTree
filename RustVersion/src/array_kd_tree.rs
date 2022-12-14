use std::{fs::File, io::Write, time::Instant};

use crate::array_particle::*;

use super::array_particle::Particle;

const MAX_PARTS: usize = 7;
const THETA: f64 = 0.3;
const NEGS: [usize; MAX_PARTS] = [usize::MAX; MAX_PARTS];

#[derive(Clone, Copy)]
pub struct KDTree {
    // For leaves
    num_parts: usize,
    particles: [usize; MAX_PARTS],

    // For internal nodes
    split_dim: usize,
    split_val: f64,
    m: f64,
    cm: [f64; 3],
    size: f64,
    left: usize,
    right: usize,
}

impl KDTree {
    pub fn leaf<'a>(num_parts: usize, particles: [usize; MAX_PARTS]) -> KDTree {
        KDTree {
            num_parts: num_parts,
            particles: particles,
            split_dim: usize::MAX,
            split_val: 0.0,
            m: 0.0,
            cm: [0.0, 0.0, 0.0],
            size: 0.0,
            left: usize::MAX,
            right: usize::MAX,
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
        nodes[cur_node].num_parts = np;
        for i in 0..np {
            nodes[cur_node].particles[i] = indices[start + i]
        }
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
        nodes[cur_node].num_parts = 0;
        nodes[cur_node].split_dim = split_dim;
        nodes[cur_node].split_val = split_val;
        nodes[cur_node].m = m;
        nodes[cur_node].cm = cm;
        nodes[cur_node].size = size;
        nodes[cur_node].left = cur_node + 1;
        nodes[cur_node].right = left + 1;

        right
    }
}

fn accel_recur(cur_node: usize, p: usize, particles: &Vec<Particle>, nodes: &Vec<KDTree>) -> [f64; 3] {
    // println!("accel {}", cur_node);
    if nodes[cur_node].num_parts > 0 {
        let mut acc = [0.0, 0.0, 0.0];
        for i in 0..(nodes[cur_node].num_parts) {
            if nodes[cur_node].particles[i] != p {
                let pp_acc = calc_pp_accel(&particles[p], &particles[nodes[cur_node].particles[i]]);
                acc[0] += pp_acc[0];
                acc[1] += pp_acc[1];
                acc[2] += pp_acc[2];
            }
        }
        acc
    } else {
        let dx = particles[p].p[0] - nodes[cur_node].cm[0];
        let dy = particles[p].p[1] - nodes[cur_node].cm[1];
        let dz = particles[p].p[2] - nodes[cur_node].cm[2];
        let dist_sqr = dx * dx + dy * dy + dz * dz;
        // println!("dist = {}, size = {}", dist, nodes[cur_node].size);
        if nodes[cur_node].size * nodes[cur_node].size < THETA * THETA * dist_sqr {
            let dist = f64::sqrt(dist_sqr);
            let magi = -nodes[cur_node].m / (dist_sqr * dist);
            [dx * magi, dy * magi, dz * magi]
        } else {
            let left_acc = accel_recur(nodes[cur_node].left, p, particles, nodes);
            let right_acc = accel_recur(nodes[cur_node].right, p, particles, nodes);
            [left_acc[0] + right_acc[0], left_acc[1] + right_acc[1], left_acc[2] + right_acc[2]]
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

    for step in 0..steps {
        // if step % 100 == 0 {
        //     let elapsed_secs = time.elapsed().as_nanos() as f64 / 1e9;
        //     println!("Step = {}, duration = {}, n = {}, nodes = {}", step, elapsed_secs, bodies.len(), tree.len());
        //     time = Instant::now();
        // }
        for i in 0..bodies.len() {
            indices[i] = i;
        }
        build_tree(&mut indices, 0, bodies.len(), bodies, 0, &mut tree);
        // if step % 10 == 0 {
        //     print_tree(step, &tree, &bodies);
        // }
        for i in 0..bodies.len() {
            acc[i] = calc_accel(i, &bodies, &tree);
        }
        for i in 0..bodies.len() {
            bodies[i].v[0] += dt * acc[i][0];
            bodies[i].v[1] += dt * acc[i][1];
            bodies[i].v[2] += dt * acc[i][2];
            let dx = dt * bodies[i].v[0];
            let dy = dt * bodies[i].v[1];
            let dz = dt * bodies[i].v[2];
            bodies[i].p[0] += dx;
            bodies[i].p[1] += dy;
            bodies[i].p[2] += dz;
            acc[i][0] = 0.0;
            acc[i][1] = 0.0;
            acc[i][2] = 0.0;
        }
    }
}

fn print_tree(step: i64, tree: &Vec<KDTree>, particles: &Vec<Particle>) -> std::io::Result<()> {
    let mut file = File::create(format!("tree{}.txt", step))?;

    file.write_fmt(format_args!("{}\n", tree.len()))?;
    for n in tree {
        if n.num_parts > 0 {
            file.write_fmt(format_args!("L {}\n", n.num_parts))?;
            for i in 0..n.num_parts {
                let p = n.particles[i];
                file.write_fmt(format_args!(
                    "{} {} {}\n",
                    particles[p].p[0], particles[p].p[1], particles[p].p[2]
                ))?;
            }
        } else {
            file.write_fmt(format_args!(
                "I {} {} {} {}\n",
                n.split_dim, n.split_val, n.left, n.right
            ))?;
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
    if nodes[node].num_parts > 0 {
        for index in 0..nodes[node].num_parts {
            let i = nodes[node].particles[index];
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
    } else {
        let split_dim = nodes[node].split_dim;
        let tmin = min[split_dim];
        let tmax = max[split_dim];
        max[split_dim] = nodes[node].split_val;
        recur_test_tree_struct(nodes[node].left, nodes, particles, min, max);
        max[split_dim] = tmax;
        min[split_dim] = nodes[node].split_val;
        recur_test_tree_struct(nodes[node].right, nodes, particles, min, max);
        min[split_dim] = tmin;
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
        assert_eq!(node_vec[0].num_parts, parts.len());
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
        assert_eq!(node_vec[0].num_parts, 0);
        assert_eq!(node_vec[1].num_parts + node_vec[2].num_parts, 12);
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
