#![feature(portable_simd)]

//mod kd_tree_box;
mod simd_particle;
mod simd_kd_tree;
mod array_particle;
mod array_kd_tree;

// use clap::Parser;
use std::{env, time::Instant};

// #[derive(Parser, Debug)]
// #[command(author, version, about, long_about = None)]
// struct Args {
//     /// Number of particules to generate.
//     #[arg(short, long)]
//     number: usize,

//     /// Number of steps to run the simulation.
//     #[arg(short, long, default_value_t = 1)]
//     steps: i64,
// }

fn main() {
    // let args = Args::parse();
    let args: Vec<String> = env::args().collect();

    let number = args[1].parse::<usize>().unwrap();
    let steps = args[2].parse::<i64>().unwrap();

    let dt = 1e-3; // * 2.0 * std::f64::consts::PI;

    let start = Instant::now();
    simd_kd_tree::simple_sim(
        &mut simd_particle::circular_orbits(number),
        dt,
        steps,
    );
    println!("{}", start.elapsed().as_nanos() as f64 / 1e9);
    let start = Instant::now();
    array_kd_tree::simple_sim(&mut array_particle::circular_orbits(number), dt, steps);
    println!("{}", start.elapsed().as_nanos() as f64 / 1e9);
}
