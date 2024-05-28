use std::time::Instant;

use rust_kdtree_nbody::quickstat::*;


fn main() {
  // Sequential
  let n = 500000000;
  let mut vals = vec!();
  let mut indices = vec!();
  for i in 0..n {
      vals.push(fastrand::f64());
      indices.push(i);
  }
  let start = 0;
  let end = n;
  let goal = fastrand::usize(start..end);
  for i in 0..n {
      vals[i] = fastrand::f64();
      indices[i] = i;
  }
  let pre = Instant::now();
  quickstat_index(&mut indices[start..end], goal - start, |i1, i2| vals[i1] < vals[i2]);
  let post = Instant::now();
  eprintln!(
    "Sequential runtime = {}",
    (post - pre).as_secs_f64());

  // Parallel
  let mut vals = vec!();
  let mut indices = vec!();
  let mut buffer: Vec<usize> = vec!();
  for i in 0..n {
      vals.push(fastrand::f64());
      indices.push(i);
      buffer.push(0);
  }
  let start = 0;
  let end = n;
  let goal = fastrand::usize(start..end);
  for i in 0..n {
      vals[i] = fastrand::f64();
      indices[i] = i;
  }
  let pre = Instant::now();
  quickstat_index_par(&mut indices[start..end], &mut buffer, goal - start, |i1, i2| vals[i1] < vals[i2]);
  let post = Instant::now();
  eprintln!(
    "Parallel runtime = {}",
    (post - pre).as_secs_f64());
}
