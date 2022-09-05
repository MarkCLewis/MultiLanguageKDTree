from __future__ import annotations

import numpy as np
from random import randint
from dataclasses import dataclass
from numpy import sqrt

from particle import f64x3, Particle, calc_pp_accel
import particle

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
    cm: f64x3
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


# Returns the index of the last Node used in the construction.
def build_tree(
    indices: list[int],
    start: int,
    end: int,
    particles: list[Particle],
    cur_node: int,
    nodes: list[KDTree],
) -> int:
    # println!("start = {} end = {} cur_node = {}", start, end, cur_node)
    np1 = end - start
    # println!("s = {}, e = {}, cn = {}", start, end, cur_node)
    if np1 <= MAX_PARTS:
        if cur_node >= len(nodes):
            diff = cur_node + 1 - len(nodes)
            nodes.extend(KDTree.leaf(0, []) for _ in range(diff))
        nodes[cur_node].num_parts = np1
        for i in range(np1):
            nodes[cur_node].particles[i] = indices[start + i]
        return cur_node
    else:
        # Pick split dim and value
        min = np.ones(3, dtype=np.float64) * 1e100
        max = np.ones(3, dtype=np.float64) * -1e100
        m = 0.0
        cm = np.zeros(3, dtype=np.float64)
        for i in range(start, end):
            m += particles[indices[i]].m
            cm += particles[indices[i]].m * particles[indices[i]].p
            min = np.minimum(min, particles[indices[i]].p)
            max = np.maximum(max, particles[indices[i]].p)

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
            # TODO: verify edge case: is this equivalent to fastrand::usize(s..e)?
            pivot = randint(s, e)
            tmp = indices[s]
            indices[s] = indices[pivot]
            indices[pivot] = tmp

            low = s + 1
            high = e - 1
            while low <= high:
                if particles[indices[low]].p[split_dim] < particles[indices[s]].p[split_dim]:
                    low += 1
                else:
                    tmp = indices[low]
                    indices[low] = indices[high]
                    indices[high] = tmp
                    high -= 1

            tmp = indices[s]
            indices[s] = indices[high]
            indices[high] = tmp

            if high < mid:
                s = high + 1
            elif high > mid:
                e = high
            else:
                s = e

        split_val = particles[indices[mid]].p[split_dim]

        # Recurse on children and build this node.
        left = build_tree(indices, start, mid, particles, cur_node + 1, nodes)
        right = build_tree(indices, mid, end, particles, left + 1, nodes)

        if cur_node >= len(nodes):
            nodes.extend(KDTree.leaf(0, [])
                         for _ in range(cur_node + 1 - len(nodes)))

        nodes[cur_node].num_parts = 0
        nodes[cur_node].split_dim = split_dim
        nodes[cur_node].split_val = split_val
        nodes[cur_node].m = m
        nodes[cur_node].cm = cm
        nodes[cur_node].size = size
        nodes[cur_node].left = cur_node+1
        nodes[cur_node].right = left+1

        return right


def accel_recur(cur_node: int, p: int, particles: list[Particle], nodes: list[KDTree]) -> f64x3:
    # println!("accel {}", cur_node)
    if nodes[cur_node].num_parts > 0:
        acc = np.zeros(3, dtype=np.float64)
        for i in range(nodes[cur_node].num_parts):
            if nodes[cur_node].particles[i] != p:
                acc += calc_pp_accel(particles[p],
                                     particles[nodes[cur_node].particles[i]])

        return acc
    else:
        dp = particles[p].p - nodes[cur_node].cm
        dist_sqr = dp @ dp
        # println!("dist = {}, size = {}", dist, nodes[cur_node].size)
        if nodes[cur_node].size * nodes[cur_node].size < THETA * THETA * dist_sqr:
            dist = sqrt(dist_sqr)
            magi = -nodes[cur_node].m / (dist_sqr*dist)
            return dp * magi
        else:
            return accel_recur(nodes[cur_node].left, p, particles, nodes) \
                + accel_recur(nodes[cur_node].right, p, particles, nodes)


def calc_accel(p: int, particles: list[Particle], nodes: list[KDTree]) -> f64x3:
    return accel_recur(0, p, particles, nodes)


def simple_sim(bodies: list[Particle], dt: float, steps: int) -> None:
    # dt_vec = f64x4:: splat(dt)
    acc = np.zeros((len(bodies), 3), dtype=np.float64)

    # time = Instant: : now()
    tree = allocate_node_vec(len(bodies))
    indices: list[int] = list(range(len(bodies)))

    for step in range(steps):
        # if step % 100 == 0 {
        # elapsed_secs = time.elapsed().as_nanos() as f64 / 1e9
        # println!("Step = {}, duration = {}, n = {}, nodes = {}", step, elapsed_secs, len(bodies), len(tree))
        #     time = Instant:: now()
        # }
        indices = list(range(len(bodies)))

        build_tree(indices, 0, len(bodies), bodies, 0, tree)
        # if step % 100 == 0 {
        # print_tree(step, tree, bodies)
        # }
        for i in range(len(bodies)):
            acc[i] = calc_accel(i, bodies, tree)

        for i in range(len(bodies)):
            bodies[i].v += dt * acc[i]
            dp = dt * bodies[i].v
            bodies[i].p += dp
            acc[i] = np.zeros(3, dtype=np.float64)


def print_tree(step: int, tree: list[KDTree], particles: list[Particle]) -> None:
    with open(f"tree{step}.txt", 'w') as file:
        file.write(f'{len(particles)}\n')
        for n in tree:
            if n.num_parts > 0:
                file.write(f"L {n.num_parts}\n")
                for i in range(n.num_parts):
                    p = n.particles[i]
                    file.write(
                        f"{particles[p].p[0]} {particles[p].p[1]} { particles[p].p[2]}\n")
            else:
                file.write(
                    f"I {n.split_dim} {n.split_val} {n.left} {n.right}\n")


def recur_test_tree_struct(
    node: int,
    nodes: list[KDTree],
    particles: list[Particle],
    min: f64x3,
    max: f64x3,
):
    if nodes[node].num_parts > 0:
        for index in range(nodes[node].num_parts):
            i = nodes[node].particles[index]
            for dim in range(2):
                assert particles[i].p[dim] >= min[dim], "Particle dim {} is below min. i={} p={} min={}".format(
                    dim,
                    i,
                    particles[i].p[dim],
                    min[dim])

                assert particles[i].p[dim] < max[dim], "Particle dim {} is above max. i={} p={} max={}".format(
                    dim,
                    i,
                    particles[i].p[dim],
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