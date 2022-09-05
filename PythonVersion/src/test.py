import numpy as np

from particle import two_bodies, circular_orbits
from kd_tree import build_tree, allocate_node_vec, recur_test_tree_struct, simple_sim


def test_single_node():
    parts = two_bodies()
    node_vec = allocate_node_vec(len(parts))
    assert len(node_vec) == 2
    indices = list(range(len(parts)))
    build_tree(indices, 0, len(parts), parts, 0, node_vec)
    assert node_vec[0].num_parts == len(parts)


def test_two_leaves():
    parts = circular_orbits(11)
    node_vec = allocate_node_vec(len(parts))
    assert len(node_vec) == 4
    indices = list(range(len(parts)))
    build_tree(indices, 0, len(parts), parts, 0, node_vec)
    recur_test_tree_struct(
        0,
        node_vec,
        parts,
        np.ones(3, dtype=np.float64) * -1e100,
        np.ones(3, dtype=np.float64) * 1e100,
    )
    assert node_vec[0].num_parts == 0
    assert node_vec[0].split_dim == 0
    assert node_vec[1].num_parts == 5
    assert node_vec[2].num_parts == 6


def test_big_solar():
    parts = circular_orbits(5000)
    node_vec = allocate_node_vec(len(parts))
    indices = list(range(len(parts)))
    build_tree(indices, 0, len(parts), parts, 0, node_vec)
    recur_test_tree_struct(
        0,
        node_vec,
        parts,
        np.ones(3, dtype=np.float64) * (-1e100),
        np.ones(3, dtype=np.float64) * (1e100),
    )


def test_big_solar_with_steps():
    parts = circular_orbits(5000)
    simple_sim(parts, 1e-3, 100)

    node_vec = allocate_node_vec(len(parts))
    indices = list(range(len(parts)))
    build_tree(indices, 0, len(parts), parts, 0, node_vec)
    recur_test_tree_struct(
        0,
        node_vec,
        parts,
        np.ones(3, dtype=np.float64) * (-1e100),
        np.ones(3, dtype=np.float64) * (1e100),
    )
