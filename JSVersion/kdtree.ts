import { closeSync, open, openSync, writeSync } from "fs";
import { calc_pp_accel, Particle, Vector3 } from "./particle";


const MAX_PARTS: number = 7
const THETA: number = 0.3;

export class KDTree {
  // For leaves
  particles: number[]

  // For internal nodes
  split_dim: number
  split_val: number
  m: number
  cm: Vector3
  size: number
  left: number
  right: number

  private constructor() {
    this.particles = []

    this.split_dim = Infinity
    this.split_val = 0
    this.m = 0
    this.cm = new Vector3(0, 0, 0)
    this.size = 0
    this.left = Infinity
    this.right = Infinity
  }

  static leaf(p: number[]) {
    let node = new KDTree()
    node.particles = p
    return node
  }
}

export function allocate_node_vec(num_parts: number): KDTree[] {
  let num_nodes = 2 * (num_parts / (MAX_PARTS - 1) + 1) | 0 // n | 0 casts n to int
  let ret = []
  for (let index = 0; index < num_nodes; index++) {
    ret.push(KDTree.leaf([]))

  }
  return ret
}

// Returns the index of the last Node used in the construction.
export function build_tree(
  indices: number[],
  start: number,
  end: number,
  particles: Particle[],
  cur_node: number,
  nodes: KDTree[],
): number {
  // println!("start = {} end = {} cur_node = {}", start, end, cur_node);
  let np = end - start;
  // println!("s = {}, e = {}, cn = {}", start, end, cur_node);
  if (np <= MAX_PARTS) {
    while (cur_node >= nodes.length) {
      nodes.push(KDTree.leaf([]))
    }
    // nodes[cur_node].num_parts = np;
    for (let i = 0; i < np; i++) {
      nodes[cur_node].particles[i] = indices[start + i]

    }

    return cur_node
  } else {
    // Pick split dim and value
    let min = new Vector3(1e100, 1e100, 1e100)
    let max = new Vector3(-1e100, -1e100, -1e100)
    let m = 0.0
    let cm = new Vector3(0, 0, 0)
    for (let i = start; i < end; i++) {
      m += particles[indices[i]].m
      cm.addUpdate(particles[indices[i]].p.scale(particles[indices[i]].m))
      min.minUpdate(particles[indices[i]].p)
      max.maxUpdate(particles[indices[i]].p)
    }
    cm.scaleUpdate(1 / m)
    let split_dim = 0;
    for (let dim = 1; dim < 3; dim++) {
      if (max.get(dim) - min.get(dim) > max.get(split_dim) - min.get(split_dim)) {
        split_dim = dim
      }
    }
    let size = max.get(split_dim) - min.get(split_dim)

    // Partition particles on split_dim
    let mid = ((start + end) / 2) | 0
    let s = start
    let e = end
    while (s + 1 < e) {
      let pivot = ((Math.random() * (e - s)) + s) | 0

      let tmp = indices[s]
      indices[s] = indices[pivot]
      indices[pivot] = tmp

      let low = s + 1
      let high = e - 1
      while (low <= high) {
        if (particles[indices[low]].p.get(split_dim) < particles[indices[s]].p.get(split_dim)) {
          low += 1;
        } else {
          let tmp = indices[low]
          indices[low] = indices[high]
          indices[high] = tmp

          high -= 1
        }
      }
      tmp = indices[s]
      indices[s] = indices[high]
      indices[high] = tmp

      if (high < mid) {
        s = high + 1;
      } else if (high > mid) {
        e = high;
      } else {
        s = e;
      }
    }
    let split_val = particles[indices[mid]].p.get(split_dim)

    // Recurse on children and build this node.
    let left = build_tree(indices, start, mid, particles, cur_node + 1, nodes);
    let right = build_tree(indices, mid, end, particles, left + 1, nodes);

    while (cur_node >= nodes.length) {
      nodes.push(KDTree.leaf([]));
    }
    nodes[cur_node].split_dim = split_dim;
    nodes[cur_node].split_val = split_val;
    nodes[cur_node].m = m;
    nodes[cur_node].cm = cm;
    nodes[cur_node].size = size;
    nodes[cur_node].left = cur_node + 1;
    nodes[cur_node].right = left + 1;

    return right
  }
}

function accel_recur(cur_node: number, p: number, particles: Particle[], nodes: KDTree[]): Vector3 {
  // println!("accel {}", cur_node);
  if (nodes[cur_node].particles.length > 0) {
    let acc = new Vector3(0, 0, 0);
    for (let i = 0; i < nodes[cur_node].particles.length; i++) {
      if (nodes[cur_node].particles[i] != p) {
        acc.addUpdate(calc_pp_accel(particles[p], particles[nodes[cur_node].particles[i]]))
      }
    }
    return acc
  } else {
    let dp = particles[p].p.sub(nodes[cur_node].cm)
    let dist_sqr = dp.magSquared()
    // println!("dist = {}, size = {}", dist, nodes[cur_node].size);
    if (nodes[cur_node].size * nodes[cur_node].size < THETA * THETA * dist_sqr) {
      let dist = Math.sqrt(dist_sqr);
      let magi = -nodes[cur_node].m / (dist_sqr * dist);
      dp.scale(magi)
      return dp
    } else {
      return accel_recur(nodes[cur_node].left, p, particles, nodes).add(
        accel_recur(nodes[cur_node].right, p, particles, nodes))
    }
  }
}

export function calc_accel(p: number, particles: Particle[], nodes: KDTree[]): Vector3 {
  return accel_recur(0, p, particles, nodes)
}

export function simple_sim(bodies: Particle[], dt: number, steps: number) {
  let acc = []
  for (let i = 0; i < bodies.length; i++) {
    acc.push(new Vector3(0, 0, 0))
  }
  // let mut time = Instant::now();
  let tree = allocate_node_vec(bodies.length)
  let indices = []
  for (let i = 0; i < bodies.length; i++) {
    indices.push(i)
  }

  for (let step = 0; step < steps; step++) {
    // if step % 100 == 0 {
    //     let elapsed_secs = time.elapsed().as_nanos() as f64 / 1e9;
    //     println!("Step = {}, duration = {}, n = {}, nodes = {}", step, elapsed_secs, bodies.len(), tree.len());
    //     time = Instant::now();
    // }
    for (let i = 0; i < bodies.length; i++) {
      indices[i] = i
    }
    build_tree(indices, 0, bodies.length, bodies, 0, tree)
    if (step % 10 === 0) {
      print_tree(step, tree, bodies)
    }
    for (let i = 0; i < bodies.length; i++) {
      acc[i] = calc_accel(i, bodies, tree)
    }
    for (let i = 0; i < bodies.length; i++) {
      bodies[i].v.addUpdate(acc[i].scale(dt))
      bodies[i].p.addUpdate(bodies[i].v.scale(dt))
      acc[i] = new Vector3(0, 0, 0)
    }
  }
}

function print_tree(step: number, tree: KDTree[], particles: Particle[]) {
  const file = openSync(`tree${step}.txt`, 'w')

  writeSync(file, `${tree.length}\n`);

  tree.forEach(n => {
    if (n.particles.length > 0) {
      writeSync(file, `L ${n.particles.length}\n`)

      n.particles.forEach(p => {
        writeSync(file,
          `${particles[p].p.x} ${particles[p].p.y} ${particles[p].p.z}\n`,
        )
      })

    } else {
      writeSync(file, `I ${n.split_dim} ${n.split_val} ${n.left} ${n.right}\n`)
    }
  })

  closeSync(file)
}

// TODO: Rust -> TS
// fn recur_test_tree_struct(
//   node: usize,
//   nodes: & Vec<KDTree>,
//   particles: & Vec<Particle>,
//   mut min: f64x4,
//   mut max: f64x4,
// ) {
//   if nodes[node].num_parts > 0 {
//     for index in 0..nodes[node].num_parts {
//       let i = nodes[node].particles[index];
//       for dim in 0..2 {
//         assert!(
//           particles[i].p[dim] >= min[dim],
//           "Particle dim {} is below min. i={} p={} min={}",
//           dim,
//           i,
//           particles[i].p[dim],
//           min[dim]
//         );
//         assert!(
//           particles[i].p[dim] < max[dim],
//           "Particle dim {} is above max. i={} p={} max={}",
//           dim,
//           i,
//           particles[i].p[dim],
//           max[dim]
//         );
//       }
//     }
//   } else {
//     let split_dim = nodes[node].split_dim;
//     let tmin = min[split_dim];
//     let tmax = max[split_dim];
//     max[split_dim] = nodes[node].split_val;
//     recur_test_tree_struct(nodes[node].left, nodes, particles, min, max);
//     max[split_dim] = tmax;
//     min[split_dim] = nodes[node].split_val;
//     recur_test_tree_struct(nodes[node].right, nodes, particles, min, max);
//     min[split_dim] = tmin;
//   }
// }

// #[cfg(test)]
// mod tests {
//     use crate:: { kd_tree, simd_particle };
//     use core_simd::*;

//     #[test]
//     fn single_node() {
//     let parts = simd_particle:: two_bodies();
//     let mut node_vec = kd_tree:: allocate_node_vec(parts.len());
//     assert_eq!(node_vec.len(), 2);
//     let mut indices: Vec<usize> = (0..parts.len()).collect();
//     kd_tree:: build_tree(& mut indices, 0, parts.len(), & parts, 0, & mut node_vec);
//     assert_eq!(node_vec[0].num_parts, parts.len());
//   }

//     #[test]
//     fn two_leaves() {
//     let parts = simd_particle:: circular_orbits(11);
//     let mut node_vec = kd_tree:: allocate_node_vec(parts.len());
//     assert_eq!(node_vec.len(), 6);
//     let mut indices: Vec<usize> = (0..parts.len()).collect();
//     kd_tree:: build_tree(& mut indices, 0, parts.len(), & parts, 0, & mut node_vec);
//     kd_tree:: recur_test_tree_struct(
//       0,
//             & node_vec,
//             & parts,
//       f64x4:: splat(-1e100),
//       f64x4:: splat(1e100),
//     );
//     assert_eq!(node_vec[0].num_parts, 0);
//     assert_eq!(node_vec[1].num_parts + node_vec[2].num_parts, 12);
//   }

//     #[test]
//     fn big_solar() {
//     let parts = simd_particle:: circular_orbits(5000);
//     let mut node_vec = kd_tree:: allocate_node_vec(parts.len());
//     let mut indices: Vec<usize> = (0..parts.len()).collect();
//     kd_tree:: build_tree(& mut indices, 0, parts.len(), & parts, 0, & mut node_vec);
//     kd_tree:: recur_test_tree_struct(
//       0,
//             & node_vec,
//             & parts,
//       f64x4:: splat(-1e100),
//       f64x4:: splat(1e100),
//     );
//   }

//     #[test]
//     fn big_solar_with_steps() {
//     let mut parts = simd_particle:: circular_orbits(5000);
//     kd_tree:: simple_sim(& mut parts, 1e-3, 10);

//     let mut node_vec = kd_tree:: allocate_node_vec(parts.len());
//     let mut indices: Vec<usize> = (0..parts.len()).collect();
//     kd_tree:: build_tree(& mut indices, 0, parts.len(), & parts, 0, & mut node_vec);
//     kd_tree:: recur_test_tree_struct(
//       0,
//             & node_vec,
//             & parts,
//       f64x4:: splat(-1e100),
//       f64x4:: splat(1e100),
//     );
//   }
// }
