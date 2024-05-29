use std::time::Instant;

use rust_kdtree_nbody::array_kd_tree::*;
use rust_kdtree_nbody::array_particle::*;

fn main() {
  let parts = circular_orbits(10000000);
  let mut node_vec = allocate_node_vec(parts.len());
  let mut indices: Vec<usize> = (0..parts.len()).collect();
  let mut buffer: Vec<usize> = (0..parts.len()).collect();
  let pre = Instant::now();
  build_tree(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
  let post = Instant::now();
  eprintln!("Sequential Runtime = {}", (post - pre).as_secs_f64());
  // let pre = Instant::now();
  // build_tree_par1(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
  // let post = Instant::now();
  // eprintln!("Par1 Runtime = {}", (post - pre).as_secs_f64());
  let pre = Instant::now();
  build_tree_par1_chunk(&mut indices, 0, parts.len(), &parts, 0, &mut node_vec);
  let post = Instant::now();
  eprintln!("Par1 Chunk Runtime = {}", (post - pre).as_secs_f64());
  let mut indices: Vec<usize> = (0..parts.len()).collect();
  let pre = Instant::now();
  build_tree_par2(&mut indices, &mut buffer, 0, parts.len(), &parts, 0, &mut node_vec);
  let post = Instant::now();
  eprintln!("Par2 Runtime = {}", (post - pre).as_secs_f64());
  let mut indices: Vec<usize> = (0..parts.len()).collect();
  let pre = Instant::now();
  build_tree_par3(&mut indices, 0, &parts, &mut node_vec, 1);
  let post = Instant::now();
  eprintln!("Par3 Runtime = {}", (post - pre).as_secs_f64());
  let pre = Instant::now();
  build_tree_par3_chunk(&mut indices, 0, &parts, &mut node_vec, 1);
  let post = Instant::now();
  eprintln!("Par3 Chunk Runtime = {}", (post - pre).as_secs_f64());
  let mut indices: Vec<usize> = (0..parts.len()).collect();
  let pre = Instant::now();
  build_tree_par4(&mut indices, 0, &parts, &mut node_vec, 1);
  let post = Instant::now();
  eprintln!("Par4 Runtime = {}", (post - pre).as_secs_f64());
}