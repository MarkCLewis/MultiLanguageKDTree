#!/usr/bin/env python3.10

from sys import argv
from time import time

import kd_tree
import particle

if __name__ == '__main__':
    n = int(argv[1])
    steps = int(argv[2])

    dt = 1e-3

    start = time()
    kd_tree.simple_sim(particle.circular_orbits(n), dt, steps)
    print(f"{time() - start}")
