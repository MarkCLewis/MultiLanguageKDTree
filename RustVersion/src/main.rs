// #![feature(portable_simd)]

// mod simd_particle;
// mod kd_tree;
mod array_particle;
mod array_kd_tree;

use std::time::Instant;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let n: usize = args[1].parse().unwrap();
    let steps: i64 = args[2].parse().unwrap();

    let dt = 1e-3; // * 2.0 * std::f64::consts::PI;

    let start = Instant::now();
    // kd_tree::simple_sim(&mut simd_particle::circular_orbits(n), dt, steps);
    array_kd_tree::simple_sim(&mut array_particle::circular_orbits(n), dt, steps);
    println!("{}", start.elapsed().as_nanos() as f64 / 1e9);
}
