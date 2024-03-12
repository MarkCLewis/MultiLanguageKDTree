use std::{fs::File, io::Write};

use crate::array_particle::*;

const MAX_PARTS: usize = 7;
const THETA: f64 = 0.3;
const NEGS: [usize; MAX_PARTS] = [usize::MAX; MAX_PARTS];

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
        left: Box<KDTree>,
        right: Box<KDTree>
    }
}


// Returns the index of the last Node used in the construction.
pub fn build_tree(
    indices: &mut Vec<usize>,
    start: usize,
    end: usize,
    particles: &Vec<Particle>,
) -> Box<KDTree> {
    // println!("start = {} end = {} cur_node = {}", start, end, cur_node);
    let np = end - start;
    // println!("s = {}, e = {}, cn = {}", start, end, cur_node);
    if np <= MAX_PARTS {
        let mut parts = [0; MAX_PARTS];
        for i in 0..np {
            parts[i] = indices[start + i]
        }
        Box::new(KDTree::Leaf { num_parts: np, leaf_parts: parts })
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
            use std::cmp::Ordering::*;
            match high.cmp(&mid) {
                Less => s = high + 1,
                Greater => e = high,
                Equal => s = e,
            }
        }
        let split_val = particles[indices[mid]].p[split_dim];

        // Recurse on children and build this node.
        let left = build_tree(indices, start, mid, particles);
        let right = build_tree(indices, mid, end, particles);

        Box::new(KDTree::Internal { split_dim, split_val, m, cm, size, left, right })
    }
}

fn calc_accel(cur_node: &Box<KDTree>, p: usize, particles: &Vec<Particle>) -> [f64; 3] {
    // println!("accel {}", cur_node);
    match &**cur_node {
        KDTree::Leaf { num_parts, leaf_parts} => {
            let mut acc = [0.0, 0.0, 0.0];
            for i in 0..*num_parts {
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
            // println!("dist = {}, size = {}", dist, size);
            if size * size < THETA * THETA * dist_sqr {
                let dist = f64::sqrt(dist_sqr);
                let magi = -m / (dist_sqr * dist);
                [dx * magi, dy * magi, dz * magi]
            } else {
                let left_acc = calc_accel(&left, p, particles);
                let right_acc = calc_accel(&right, p, particles);
                [left_acc[0] + right_acc[0], left_acc[1] + right_acc[1], left_acc[2] + right_acc[2]]
            }
        }
    }
}

pub fn simple_sim(bodies: &mut Vec<Particle>, dt: f64, steps: i64) {
    let mut acc = Vec::new();
    for _ in 0..bodies.len() {
        acc.push([0.0, 0.0, 0.0])
    }
    // let mut time = Instant::now();
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

        let root = build_tree(&mut indices, 0, bodies.len(), bodies);
        // if step % 10 == 0 {
        //     print_tree(step, &tree, &bodies);
        // }
        for (index, item) in acc.iter_mut().enumerate().take(bodies.len()) {
            *item = calc_accel(&root, index, bodies)
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

#[allow(dead_code)]
fn print_tree(step: i64, tree: &Box<KDTree>, particles: &Vec<Particle>) -> std::io::Result<()> {
    let mut file = File::create(format!("tree{}.txt", step))?;

    // file.write_fmt(format_args!("{}\n", tree.len()))?;
    // for n in tree {
    //     match n {
    //         KDTree::Leaf { num_parts, leaf_parts } => {
    //             file.write_fmt(format_args!("L {}\n", num_parts))?;
    //             for i in 0..*num_parts {
    //                 let p = leaf_parts[i];
    //                 file.write_fmt(format_args!(
    //                     "{} {} {}\n",
    //                     particles[p].p[0], particles[p].p[1], particles[p].p[2]
    //                 ))?;
    //             }
    //         }
    //         KDTree::Internal { split_dim, split_val, left, right, .. } => {
    //             file.write_fmt(format_args!(
    //                 "I {} {} {} {}\n",
    //                 split_dim, split_val, left, right
    //             ))?;
    //         }
    //     }
    // }

    Ok(())
}

fn recur_test_tree_struct(
    node: &Box<KDTree>,
    particles: &Vec<Particle>,
    mut min: [f64; 3],
    mut max: [f64; 3],
) {
    match &**node {
        KDTree::Leaf { num_parts, leaf_parts } => {
            for index in 0..*num_parts {
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
            let tmin = min[*split_dim];
            let tmax = max[*split_dim];
            max[*split_dim] = *split_val;
            recur_test_tree_struct(left, particles, min, max);
            max[*split_dim] = tmax;
            min[*split_dim] = *split_val;
            recur_test_tree_struct(right, particles, min, max);
            min[*split_dim] = tmin;
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