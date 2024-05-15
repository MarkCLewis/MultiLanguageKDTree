#include<iostream>
#include<omp.h>

#include "particle.h"
#include "kdtree.h"

int main(int argc, char *argv[]) {
	if (argc < 3) {
		std::cout << "Specify a number of particles and number of steps.\n";
		return 1;
	}
	std::cout << "Running sim.\n";

	int steps = atoi(argv[1]);
	int n = atoi(argv[2]);
        omp_set_num_threads(std::atoi(argv[3]));

	double dt = 1e-3; // * 2.0 * std::f64::consts::PI;

	// let start = Instant::now();
	auto particles = circular_orbits(n);
	simple_sim(particles, dt, steps);
	// println!("{}", start.elapsed().as_nanos() as f64 / 1e9);

	return 0;
}
