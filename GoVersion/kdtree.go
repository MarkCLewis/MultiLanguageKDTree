package main

import (
	"math"
	"math/rand"
)

const (
	MAX_PARTS int     = 7
	THETA     float64 = 0.3
)

type KDTree struct {
	// For leaves
	num_parts int
	particles [MAX_PARTS]int

	// For internal nodes
	split_dim int
	split_val float64
	m         float64
	cm        [3]float64
	size      float64
	left      int
	right     int
}

func empty_leaf() KDTree {
	return KDTree{
		0,
		[MAX_PARTS]int{},
		-1,
		0.0,
		0.0,
		[3]float64{0.0, 0.0, 0.0},
		0.0,
		-1,
		-1,
	}
}

func allocate_node_vec(num_parts int) []KDTree {
	num_nodes := 2 * (num_parts/(MAX_PARTS-1) + 1)
	ret := make([]KDTree, num_nodes)
	return ret
}

// Returns the index of the last Node used in the construction.
func build_tree(
	indices []int,
	start int,
	end int,
	particles []Particle,
	cur_node int,
	nodes []KDTree,
) int {
	// println!("start = {} end = {} cur_node = {}", start, end, cur_node);
	np := end - start
	// println!("s = {}, e = {}, cn = {}", start, end, cur_node);
	if np <= MAX_PARTS {
		for cur_node >= len(nodes) {
			nodes = append(nodes, empty_leaf())
		}
		nodes[cur_node].num_parts = np
		for i := 0; i < np; i++ {
			nodes[cur_node].particles[i] = indices[start+i]
		}
		return cur_node
	} else {
		// Pick split dim and value
		min := [3]float64{1e100, 1e100, 1e100}
		max := [3]float64{-1e100, -1e100, -1e100}
		m := 0.0
		cm := [3]float64{0.0, 0.0, 0.0}
		for i := start; i < end; i++ {
			m += particles[indices[i]].m
			cm[0] += particles[indices[i]].m * particles[indices[i]].p[0]
			cm[1] += particles[indices[i]].m * particles[indices[i]].p[1]
			cm[2] += particles[indices[i]].m * particles[indices[i]].p[2]
			min[0] = math.Min(min[0], particles[indices[i]].p[0])
			min[1] = math.Min(min[1], particles[indices[i]].p[1])
			min[2] = math.Min(min[2], particles[indices[i]].p[2])
			max[0] = math.Max(max[0], particles[indices[i]].p[0])
			max[1] = math.Max(max[1], particles[indices[i]].p[1])
			max[2] = math.Max(max[2], particles[indices[i]].p[2])
		}
		cm[0] /= m
		cm[1] /= m
		cm[2] /= m
		split_dim := 0
		for dim := 1; dim < 3; dim++ {
			if max[dim]-min[dim] > max[split_dim]-min[split_dim] {
				split_dim = dim
			}
		}
		size := max[split_dim] - min[split_dim]

		// Partition particles on split_dim
		mid := (start + end) / 2
		s := start
		e := end
		for s+1 < e {
			pivot := s + rand.Int()%(e-s)
			tmp := indices[s]
			indices[s] = indices[pivot]
			indices[pivot] = tmp
			low := s + 1
			high := e - 1
			for low <= high {
				if particles[indices[low]].p[split_dim] < particles[indices[s]].p[split_dim] {
					low += 1
				} else {
					tmp = indices[low]
					indices[low] = indices[high]
					indices[high] = tmp
					high -= 1
				}
			}
			tmp = indices[s]
			indices[s] = indices[high]
			indices[high] = tmp
			if high < mid {
				s = high + 1
			} else if high > mid {
				e = high
			} else {
				s = e
			}
		}
		split_val := particles[indices[mid]].p[split_dim]

		// Recurse on children and build this node.
		left := build_tree(indices, start, mid, particles, cur_node+1, nodes)
		right := build_tree(indices, mid, end, particles, left+1, nodes)

		for cur_node >= len(nodes) {
			nodes = append(nodes, empty_leaf())
		}
		nodes[cur_node].num_parts = 0
		nodes[cur_node].split_dim = split_dim
		nodes[cur_node].split_val = split_val
		nodes[cur_node].m = m
		nodes[cur_node].cm = cm
		nodes[cur_node].size = size
		nodes[cur_node].left = cur_node + 1
		nodes[cur_node].right = left + 1

		return right
	}
}

func accel_recur(cur_node int, p int, particles []Particle, nodes []KDTree) [3]float64 {
	// println!("accel {}", cur_node);
	if nodes[cur_node].num_parts > 0 {
		acc := [3]float64{0.0, 0.0, 0.0}
		for i := 0; i < nodes[cur_node].num_parts; i++ {
			if nodes[cur_node].particles[i] != p {
				pp_acc := Calc_pp_accel(&particles[p], &particles[nodes[cur_node].particles[i]])
				acc[0] += pp_acc[0]
				acc[1] += pp_acc[1]
				acc[2] += pp_acc[2]
			}
		}
		return acc
	} else {
		dx := particles[p].p[0] - nodes[cur_node].cm[0]
		dy := particles[p].p[1] - nodes[cur_node].cm[1]
		dz := particles[p].p[2] - nodes[cur_node].cm[2]
		dist_sqr := dx*dx + dy*dy + dz*dz
		// println!("dist = {}, size = {}", dist, nodes[cur_node].size);
		if nodes[cur_node].size*nodes[cur_node].size < THETA*THETA*dist_sqr {
			dist := math.Sqrt(dist_sqr)
			magi := -nodes[cur_node].m / (dist_sqr * dist)
			return [3]float64{dx * magi, dy * magi, dz * magi}
		} else {
			left_acc := accel_recur(nodes[cur_node].left, p, particles, nodes)
			right_acc := accel_recur(nodes[cur_node].right, p, particles, nodes)
			return [3]float64{left_acc[0] + right_acc[0], left_acc[1] + right_acc[1], left_acc[2] + right_acc[2]}
		}
	}
}

func calc_accel(p int, particles []Particle, nodes []KDTree) [3]float64 {
	return accel_recur(0, p, particles, nodes)
}

func Simple_sim(bodies []Particle, dt float64, steps int) {
	acc := make([][3]float64, len(bodies))
	tree := allocate_node_vec(len(bodies))
	indices := make([]int, len(bodies))

	for step := 0; step < steps; step++ {
		// if step % 100 == 0 {
		//     let elapsed_secs = time.elapsed().as_nanos() as f64 / 1e9;
		//     println!("Step = {}, duration = {}, n = {}, nodes = {}", step, elapsed_secs, bodies.len(), tree.len());
		//     time = Instant::now();
		// }
		for i := 0; i < len(bodies); i++ {
			indices[i] = i
		}
		build_tree(indices, 0, len(bodies), bodies, 0, tree)
		// if step % 10 == 0 {
		//     print_tree(step, &tree, &bodies);
		// }
		for i := 0; i < len(bodies); i++ {
			acc[i] = calc_accel(i, bodies, tree)
		}
		for i := 0; i < len(bodies); i++ {
			bodies[i].v[0] += dt * acc[i][0]
			bodies[i].v[1] += dt * acc[i][1]
			bodies[i].v[2] += dt * acc[i][2]
			dx := dt * bodies[i].v[0]
			dy := dt * bodies[i].v[1]
			dz := dt * bodies[i].v[2]
			bodies[i].p[0] += dx
			bodies[i].p[1] += dy
			bodies[i].p[2] += dz
			acc[i][0] = 0.0
			acc[i][1] = 0.0
			acc[i][2] = 0.0
		}
	}
}

// func print_tree(step int, tree []KDTree, particles []Particle) {
//     let mut file = File::create(format!("tree{}.txt", step))?;

//     file.write_fmt(format_args!("{}\n", tree.len()))?;
//     for n in tree {
//         if n.num_parts > 0 {
//             file.write_fmt(format_args!("L {}\n", n.num_parts))?;
//             for i in 0..n.num_parts {
//                 let p = n.particles[i];
//                 file.write_fmt(format_args!(
//                     "{} {} {}\n",
//                     particles[p].p[0], particles[p].p[1], particles[p].p[2]
//                 ))?;
//             }
//         } else {
//             file.write_fmt(format_args!(
//                 "I {} {} {} {}\n",
//                 n.split_dim, n.split_val, n.left, n.right
//             ))?;
//         }
//     }

//     Ok(())
// }

// fn recur_test_tree_struct(
//     node: usize,
//     nodes: &Vec<KDTree>,
//     particles: &Vec<Particle>,
//     mut min: [f64; 3],
//     mut max: [f64; 3],
// ) {
//     if nodes[node].num_parts > 0 {
//         for index in 0..nodes[node].num_parts {
//             let i = nodes[node].particles[index];
//             for dim in 0..2 {
//                 assert!(
//                     particles[i].p[dim] >= min[dim],
//                     "Particle dim {} is below min. i={} p={} min={}",
//                     dim,
//                     i,
//                     particles[i].p[dim],
//                     min[dim]
//                 );
//                 assert!(
//                     particles[i].p[dim] < max[dim],
//                     "Particle dim {} is above max. i={} p={} max={}",
//                     dim,
//                     i,
//                     particles[i].p[dim],
//                     max[dim]
//                 );
//             }
//         }
//     } else {
//         let split_dim = nodes[node].split_dim;
//         let tmin = min[split_dim];
//         let tmax = max[split_dim];
//         max[split_dim] = nodes[node].split_val;
//         recur_test_tree_struct(nodes[node].left, nodes, particles, min, max);
//         max[split_dim] = tmax;
//         min[split_dim] = nodes[node].split_val;
//         recur_test_tree_struct(nodes[node].right, nodes, particles, min, max);
//         min[split_dim] = tmin;
//     }
// }

// #[cfg(test)]
// mod tests {
//     use crate::{array_kd_tree, array_particle};

//     #[test]
//     fn single_node() {
//         let parts = array_particle::two_bodies();
//         let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
//         assert_eq!(node_vec.len(), 2);
//         let mut indices: Vec<usize> = (0..parts.len()).collect();
//         array_kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
//         assert_eq!(node_vec[0].num_parts, parts.len());
//     }

//     #[test]
//     fn two_leaves() {
//         let parts = array_particle::circular_orbits(11);
//         let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
//         assert_eq!(node_vec.len(), 6);
//         let mut indices: Vec<usize> = (0..parts.len()).collect();
//         array_kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
//         array_kd_tree::recur_test_tree_struct(
//             0,
//             &node_vec,
//             &parts,
//             [-1e100, -1e100, -1e100],
//             [1e100, 1e100, 1e100],
//         );
//         assert_eq!(node_vec[0].num_parts, 0);
//         assert_eq!(node_vec[1].num_parts + node_vec[2].num_parts, 12);
//     }

//     #[test]
//     fn big_solar() {
//         let parts = array_particle::circular_orbits(5000);
//         let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
//         let mut indices: Vec<usize> = (0..parts.len()).collect();
//         array_kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
//         array_kd_tree::recur_test_tree_struct(
//             0,
//             &node_vec,
//             &parts,
//             [-1e100, -1e100, -1e100],
//             [1e100, 1e100, 1e100],
//         );
//     }

//     #[test]
//     fn big_solar_with_steps() {
//         let mut parts = array_particle::circular_orbits(5000);
//         array_kd_tree::simple_sim(&mut parts, 1e-3, 10);

//         let mut node_vec = array_kd_tree::allocate_node_vec(parts.len());
//         let mut indices: Vec<usize> = (0..parts.len()).collect();
//         array_kd_tree::build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
//         array_kd_tree::recur_test_tree_struct(
//             0,
//             &node_vec,
//             &parts,
//             [-1e100, -1e100, -1e100],
//             [1e100, 1e100, 1e100],
//         );
//     }
// }
