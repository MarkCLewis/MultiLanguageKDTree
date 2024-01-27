import { allocate_node_vec, build_tree, KDTree, simple_sim } from "./kdtree";
import { circular_orbits, Particle, two_bodies, Vector3 } from "./particle";
import assert from "assert";


function recur_test_tree_struct(node: number, nodes: KDTree[], particles: Particle[], min: Vector3, max: Vector3) {
  if (nodes[node].particles.length > 0) {
    for (let index = 0; index < nodes[node].particles.length; index++) {
      let i = nodes[node].particles[index];
      for (let dim = 0; dim < 2; dim++) {
        assert(particles[i].p.get(dim) >= min.get(dim))
        //   "Particle dim {} is below min. i={} p={} min={}",
        //   dim,
        //   i,
        //   particles[i].p[dim],
        //   min[dim]
        // );
        assert(particles[i].p.get(dim) < max.get(dim))
        //   "Particle dim {} is above max. i={} p={} max={}",
        //   dim,
        //   i,
        //   particles[i].p[dim],
        //   max[dim]
        // );
      }
    }
  } else {
    let split_dim = nodes[node].split_dim
    let tmin = min.get(split_dim)
    let tmax = max.get(split_dim)
    max.set(split_dim, nodes[node].split_val)
    recur_test_tree_struct(nodes[node].left, nodes, particles, min, max);
    max.set(split_dim, tmax)
    min.set(split_dim, nodes[node].split_val)
    recur_test_tree_struct(nodes[node].right, nodes, particles, min, max);
    min.set(split_dim, tmin)
  }
}

export function single_node() {
  let parts = two_bodies()
  let node_vec = allocate_node_vec(parts.length)
  assert.strictEqual(node_vec.length, 2)

  let indices: number[] = []
  for (let i = 0; i < parts.length; i++) {
    indices[i] = i;
  }

  build_tree(indices, 0, parts.length, parts, 0, node_vec)
  assert.strictEqual(node_vec[0].particles.length, parts.length)
}


export function two_leaves() {
  let parts = circular_orbits(11);
  let node_vec = allocate_node_vec(parts.length);
  assert.strictEqual(node_vec.length, 6);
  let indices = [];
  for (let i = 0; i < parts.length; i++) {
    indices.push(i);
  }
  build_tree(indices, 0, parts.length, parts, 0, node_vec);
  recur_test_tree_struct(
    0,
    node_vec,
    parts,
    new Vector3(-1e100, -1e100, -1e100),
    new Vector3(1e100, 1e100, 1e100)
  );
  assert.strictEqual(node_vec[0].particles.length, 0);
  assert.strictEqual(node_vec[1].particles.length + node_vec[2].particles.length, 12);
}

export function big_solar() {
  let parts = circular_orbits(5000);
  let node_vec = allocate_node_vec(parts.length);
  let indices = [];
  for (let i = 0; i < parts.length; i++) {
    indices.push(i);
  }
  build_tree(indices, 0, parts.length, parts, 0, node_vec);
  recur_test_tree_struct(
    0,
    node_vec,
    parts,
    new Vector3(-1e100, -1e100, -1e100),
    new Vector3(1e100, 1e100, 1e100)
  );
}

export function big_solar_with_steps() {
  let parts = circular_orbits(5000);
  simple_sim(parts, 1e-3, 10);

  let node_vec = allocate_node_vec(parts.length);
  let indices = [];
  for (let i = 0; i < parts.length; i++) {
    indices.push(i);
  }
  build_tree(indices, 0, parts.length, parts, 0, node_vec);
  recur_test_tree_struct(
    0,
    node_vec,
    parts,
    new Vector3(-1e100, -1e100, -1e100),
    new Vector3(1e100, 1e100, 1e100)
  );
}

