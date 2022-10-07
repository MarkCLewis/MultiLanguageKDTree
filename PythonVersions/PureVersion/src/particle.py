'''
Written for python3.10
'''
from __future__ import annotations

from dataclasses import dataclass
from random import random
from math import sqrt, cos, sin
# import numpy as np
# from numpy import sqrt, cos, sin, sum
# import numpy.typing as npt


# f64x3 = npt.NDArray[np.float64]


class F64x3:
    __slots__ = ['x', 'y', 'z']

    def __init__(self, x, y, z) -> None:
        self.x = x
        self.y = y
        self.z = z

    def min(self, other):
        return F64x3(min(self.x, other.x), min(self.y, other.y), min(self.z, other.z))

    def max(self, other):
        return F64x3(max(self.x, other.x), max(self.y, other.y), max(self.z, other.z))

    def __add__(self, other):
        return F64x3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return F64x3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __imul__(self, other):
        self.x *= other
        self.y *= other
        self.z *= other
        return self

    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def __mul__(self, other):
        return F64x3(self.x * other, self.y * other, self.z * other)

    def __rmul__(self, other):
        # multiplication is commutative
        return self.__mul__(other)

    def __truediv__(self, other):
        return F64x3(self.x / other, self.y / other, self.z / other)

    def __itruediv__(self, other):
        self.x /= other
        self.y /= other
        self.z /= other
        return self

    def __matmul__(self, other):
        return self.x * other.x + self.y * other.y + self.z * other.z

    def __getitem__(self, idx):
        if idx == 0:
            return self.x
        elif idx == 1:
            return self.y
        elif idx == 2:
            return self.z
        raise IndexError(f"Got index {idx}, but only 0-2 are valid")


f64x3 = F64x3


@dataclass(slots=True, init=False)
class Particle:
    def __init__(self, p: f64x3 | tuple[float, float, float], v: f64x3 | tuple[float, float, float], r: float, m: float) -> None:
        self.r = r
        self.m = m
        self.p = p if isinstance(p, F64x3) else F64x3(*p)
        self.v = v if isinstance(v, F64x3) else F64x3(*v)
        # self.p = p if isinstance(
        #     p, np.ndarray) else np.array(p, dtype=np.float64)
        # self.v = v if isinstance(
        #     v, np.ndarray) else np.array(v, dtype=np.float64)

    p: f64x3
    v: f64x3
    r: float
    m: float


def two_bodies() -> list[Particle]:
    return [
        Particle((0, 0, 0), (0, 0, 0), 1, 1),
        Particle((1, 0, 0), (0, 1, 0), 1e-4, 1e-20)
    ]


def circular_orbits(n: int) -> list[Particle]:
    particles = [Particle((0, 0, 0), (0, 0, 0), 0.00465047, 1)]

    for i in range(n):
        d = 0.1 + (i * 5 / n)
        v = sqrt(1.0 / d)
        theta = random() * 6.28
        x = d * cos(theta)
        y = d * sin(theta)
        vx = -v * sin(theta)
        vy = v * cos(theta)
        particles.append(Particle(
            (x, y, 0),
            (vx, vy, 0),
            1e-14,
            1e-7,
        ))
    return particles


def calc_pp_accel(pi: Particle, pj: Particle) -> f64x3:
    dp = pi.p - pj.p
    dp2 = dp @ dp
    dist = sqrt(dp2)
    magi = -pj.m / (dist * dist * dist)
    return dp * magi
