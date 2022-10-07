'''
Written for python3.10
'''
from __future__ import annotations

from dataclasses import dataclass
from random import random
import numpy as np
from numpy import sqrt, cos, sin, sum
import numpy.typing as npt


f64Arr = npt.NDArray[np.float64]

Particles = tuple[f64Arr, f64Arr, f64Arr, f64Arr]


def make_particle_arrays(particles: list[tuple[tuple[float, float, float], tuple[float, float, float], float, float]]) -> Particles:
    length = len(particles)
    p = np.zeros((length, 3), dtype=np.float64)
    v = np.zeros((length, 3), dtype=np.float64)
    r = np.zeros(length, dtype=np.float64)
    m = np.zeros(length, dtype=np.float64)

    for i, (p_i, v_i, r_i, m_i) in enumerate(particles):
        p[i][:] = p_i
        v[i][:] = v_i
        r[i] = r_i
        m[i] = m_i

    return p, v, r, m


def two_bodies() -> Particles:
    return make_particle_arrays([
        ((0, 0, 0), (0, 0, 0), 1, 1),
        ((1, 0, 0), (0, 1, 0), 1e-4, 1e-20)
    ])


def circular_orbits(n: int) -> Particles:
    particles = [((0, 0, 0), (0, 0, 0), 0.00465047, 1)]

    for i in range(n):
        d = 0.1 + (i * 5 / n)
        v = sqrt(1.0 / d)
        theta = random() * 6.28
        x = d * cos(theta)
        y = d * sin(theta)
        vx = -v * sin(theta)
        vy = v * cos(theta)
        particles.append((
            (x, y, 0),
            (vx, vy, 0),
            1e-14,
            1e-7,
        ))
    return make_particle_arrays(particles)


def calc_pp_accel(pi, pj, mj):
    dp = pi - pj
    dp2 = dp @ dp
    dist = sqrt(dp2)
    magi = -mj / (dist * dist * dist)
    return dp * magi
