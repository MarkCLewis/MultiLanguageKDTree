#![feature(portable_simd)]

//mod kd_tree_box;
// mod simd_particle;
mod array_particle;
mod array_kd_tree;

use clap::Parser;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Number of particules to generate.
    #[arg(short, long)]
    number: usize,

    /// Number of steps to run the simulation.
    #[arg(short, long, default_value_t = 1)]
    steps: i64,
}

fn main() {
    let args = Args::parse();

    let dt = 1e-3; // * 2.0 * std::f64::consts::PI;

    let start = Instant::now();
    // kd_tree::simple_sim(
    //     &mut simd_particle::circular_orbits(args.number),
    //     dt,
    //     args.steps,
    // );
    array_kd_tree::simple_sim(&mut array_particle::circular_orbits(args.number), dt, args.steps);
    println!("{}", start.elapsed().as_nanos() as f64 / 1e9);
}
