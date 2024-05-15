#include "kdtree.h"
#include "particle.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char *argv[]) {
  if (argc < 3) {
    printf("Specify a number of particles and number of steps.\n");
    return 1;
  }
  printf("Running sim.\n");

  int steps = atoi(argv[1]);
  int n = atoi(argv[2]);
  omp_set_num_threads(atoi(argv[3]));

  double dt = 1e-3; // * 2.0 * std::f64::consts::PI;

  // let start = Instant::now();
  Particle_array_t particles = circular_orbits(n);
  simple_sim(&particles, dt, steps);
  // println!("{}", start.elapsed().as_nanos() as f64 / 1e9);

  FREE_ARRAY(particles);

  return 0;
}
