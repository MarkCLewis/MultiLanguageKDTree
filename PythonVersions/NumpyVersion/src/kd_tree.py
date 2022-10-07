from __future__ import annotations

import numpy as np
from random import randrange
from dataclasses import dataclass
from numpy import sqrt
import numpy.typing as npt

from particle import calc_pp_accel, Particles


MAX_PARTS: int = 7
THETA: float = 0.3


@dataclass(slots=True)
class KDTree:
    # For leaves
    num_parts: int
    particles: list[int]

    # For internal nodes
    split_dim: int
    split_val: float
    m: float
    cm: tuple[float, float, float]  # de facto np.array
    size: float
    left: int
    right: int

    @staticmethod
    def leaf(num_parts: int, particles: list[int]) -> KDTree:
        return KDTree(
            num_parts,
            particles,
            None,  # type: ignore
            0.0,
            0.0,
            np.zeros(3, dtype=np.float64),
            0.0,
            None,  # type: ignore
            None  # type: ignore
        )


def allocate_node_vec(num_parts: int) -> list[KDTree]:
    num_nodes = 2 * (num_parts // (MAX_PARTS-1) + 1)
    return [KDTree.leaf(0, []) for _ in range(num_nodes)]


class System:
    indices: list[int]
    nodes: list[KDTree]

    __slots__ = ['indices', 'nodes']

    def __init__(self, n: int) -> None:
        self.indices = None
        self.nodes = allocate_node_vec(n)

    # Returns the index of the last Node used in the construction.
    def build_tree(self,
                   start: int,
                   end: int,
                   particles: Particles,
                   cur_node: int,
                   ) -> int:
        ps, vs, rs, ms = particles
        del particles

        # println!("start = {} end = {} cur_node = {}", start, end, cur_node)
        np1 = end - start
        # println!("s = {}, e = {}, cn = {}", start, end, cur_node)
        if np1 <= MAX_PARTS:
            if cur_node >= len(self.nodes):
                diff = cur_node + 1 - len(self.nodes)
                self.nodes.extend(KDTree.leaf(0, []) for _ in range(diff))
            self.nodes[cur_node].num_parts = np1

            self.nodes[cur_node].particles = [self.indices[start + i]
                                              for i in range(np1)]
            return cur_node
        else:
            # Pick split dim and value
            min = np.ones(3, dtype=np.float64) * 1e100
            max = np.ones(3, dtype=np.float64) * -1e100
            m = 0.0
            cm = np.zeros(3, dtype=np.float64)
            for i in range(start, end):
                m += ms[self.indices[i]]
                cm += ms[self.indices[i]] * \
                    ps[self.indices[i]]
                min = np.minimum(min, ps[self.indices[i]])
                max = np.maximum(max, ps[self.indices[i]])

            cm /= m
            split_dim = 0
            for dim in range(1, 2):  # FIXME: what is this
                if max[dim] - min[dim] > max[split_dim] - min[split_dim]:
                    split_dim = dim
            size = max[split_dim] - min[split_dim]

            # Partition particles on split_dim
            mid = (start + end) // 2
            s = start
            e = end
            while s + 1 < e:
                pivot = randrange(s, e)
                tmp = self.indices[s]
                try:
                    self.indices[s] = self.indices[pivot]
                except IndexError:
                    print(s, pivot, len(self.indices))
                self.indices[pivot] = tmp

                low = s + 1
                high = e - 1
                while low <= high:
                    if ps[self.indices[low]][split_dim] < ps[self.indices[s]][split_dim]:
                        low += 1
                    else:
                        tmp = self.indices[low]
                        self.indices[low] = self.indices[high]
                        self.indices[high] = tmp
                        high -= 1

                tmp = self.indices[s]
                self.indices[s] = self.indices[high]
                self.indices[high] = tmp

                if high < mid:
                    s = high + 1
                elif high > mid:
                    e = high
                else:
                    s = e

            split_val = ps[self.indices[mid]][split_dim]

            # Recurse on children and build this node.
            left = self.build_tree(start, mid,
                                   (ps, vs, rs, ms), cur_node + 1)
            right = self.build_tree(mid, end,
                                    (ps, vs, rs, ms), left + 1)

            if cur_node >= len(self.nodes):
                self.nodes.extend(KDTree.leaf(0, [])
                                  for _ in range(cur_node + 1 - len(self.nodes)))

            self.nodes[cur_node].num_parts = 0
            self.nodes[cur_node].split_dim = split_dim
            self.nodes[cur_node].split_val = split_val
            self.nodes[cur_node].m = m
            self.nodes[cur_node].cm = cm
            self.nodes[cur_node].size = size
            self.nodes[cur_node].left = cur_node+1
            self.nodes[cur_node].right = left+1

            return right


def accel_recur(cur_node: int, p: int, particles: Particles, nodes: list[KDTree]) -> npt.NDArray[np.float64]:
    ps, vs, rs, ms = particles
    del particles

    # println!("accel {}", cur_node)
    if nodes[cur_node].num_parts > 0:
        acc = np.zeros(3, dtype=np.float64)
        # print(nodes[cur_node].particles)

        # acc = calc_pp_accel(particles[p],
        #                     particles[nodes[cur_node].particles])
        for i in range(nodes[cur_node].num_parts):
            if nodes[cur_node].particles[i] != p:
                idx = nodes[cur_node].particles[i]
                acc += calc_pp_accel(ps[p],
                                     ps[idx], ms[idx])

        return acc
    else:
        dp = ps[p] - nodes[cur_node].cm
        dist_sqr = dp @ dp
        # println!("dist = {}, size = {}", dist, nodes[cur_node].size)
        if nodes[cur_node].size * nodes[cur_node].size < THETA * THETA * dist_sqr:
            dist = sqrt(dist_sqr)
            magi = -nodes[cur_node].m / (dist_sqr*dist)
            return dp * magi
        else:
            return accel_recur(nodes[cur_node].left, p, (ps, vs, rs, ms), nodes) \
                + accel_recur(nodes[cur_node].right, p,
                              (ps, vs, rs, ms), nodes)


def calc_accel(p: int, particles: Particles, nodes: list[KDTree]):
    return accel_recur(0, p, particles, nodes)


def simple_sim(bodies: Particles, dt: float, steps: int, print_steps: bool = False) -> None:
    p, v, r, m = bodies
    del bodies

    acc = np.zeros((len(p), 3), dtype=np.float64)

    sys = System(len(p))

    dv = np.zeros((len(p), 3), dtype=np.float64)

    for step in range(steps):
        if print_steps:
            print(step)
        sys.indices = [i for i in range(len(p))]

        sys.build_tree(0, len(p), (p, v, r, m), 0)
        if step % 10 == 0:
            print_tree(step, sys.nodes, (p, v, r, m))

        for i in range(len(p)):
            acc[i] = calc_accel(i, (p, v, r, m), sys.nodes)

        acc[:] *= dt  # dt * acc

        v[:] += acc
        np.multiply(dt, v, out=dv)
        p[:] += dv  # dt * v


def print_tree(step: int, tree: list[KDTree], particles: Particles) -> None:
    ps, vs, rs, ms = particles

    with open(f"tree{step}.txt", 'w') as file:
        file.write(f'{len(ps)}\n')
        for n in tree:
            if n.num_parts > 0:
                file.write(f"L {n.num_parts}\n")
                for i in range(n.num_parts):
                    p = n.particles[i]
                    file.write(
                        f"{ps[p][0]} {ps[p][1]} {ps[p][2]}\n")
            else:
                file.write(
                    f"I {n.split_dim} {n.split_val} {n.left} {n.right}\n")


def recur_test_tree_struct(
    node: int,
    nodes: list[KDTree],
    particles: Particles,
    min: tuple[float, float, float],
    max: tuple[float, float, float],
):
    ps, vs, rs, ms = particles

    if nodes[node].num_parts > 0:
        for index in range(nodes[node].num_parts):
            i = nodes[node].particles[index]
            for dim in range(2):
                assert ps[i][dim] >= min[dim], "Particle dim {} is below min. i={} p={} min={}".format(
                    dim,
                    i,
                    ps[i][dim],
                    min[dim])

                assert ps[i][dim] < max[dim], "Particle dim {} is above max. i={} p={} max={}".format(
                    dim,
                    i,
                    ps[i][dim],
                    max[dim])

    else:
        split_dim = nodes[node].split_dim
        tmin = min[split_dim]
        tmax = max[split_dim]
        max[split_dim] = nodes[node].split_val
        recur_test_tree_struct(nodes[node].left, nodes, particles, min, max)
        max[split_dim] = tmax
        min[split_dim] = nodes[node].split_val
        recur_test_tree_struct(nodes[node].right, nodes, particles, min, max)
        min[split_dim] = tmin
