import numpy as np

from particle import two_bodies, circular_orbits
from kd_tree import System, recur_test_tree_struct, simple_sim, print_tree


def test_single_node():
    parts = two_bodies()
    sys = System.from_amount(len(parts))
    assert len(sys.nodes) == 2
    sys.build_tree(0, len(parts), parts, 0)
    assert sys.nodes[0].num_parts == len(parts)


def test_two_leaves():
    parts = circular_orbits(11)
    sys = System.from_amount(len(parts))
    assert len(sys.nodes) == 4  # In the rust version, this fails also
    sys.build_tree(0, len(parts), parts, 0)
    recur_test_tree_struct(
        0,
        sys.nodes,
        parts,
        np.ones(3, dtype=np.float64) * -1e100,
        np.ones(3, dtype=np.float64) * 1e100,
    )
    assert sys.nodes[0].num_parts == 0
    assert sys.nodes[0].split_dim == 0
    assert sys.nodes[1].num_parts == 5
    assert sys.nodes[2].num_parts == 6


def test_big_solar():
    parts = circular_orbits(5000)
    sys = System.from_amount(len(parts))
    sys.build_tree(0, len(parts), parts, 0)
    recur_test_tree_struct(
        0,
        sys.nodes,
        parts,
        np.ones(3, dtype=np.float64) * (-1e100),
        np.ones(3, dtype=np.float64) * (1e100),
    )

    # # TODO: remove again
    # print_tree(0, sys.nodes, parts)


def test_big_solar_with_steps():
    parts = circular_orbits(5000)
    simple_sim(parts, 1e-3, 10, print_steps=True)

    sys = System.from_amount(len(parts))

    sys.build_tree(0, len(parts), parts, 0)
    recur_test_tree_struct(
        0,
        sys.nodes,
        parts,
        np.ones(3, dtype=np.float64) * (-1e100),
        np.ones(3, dtype=np.float64) * (1e100),
    )


if __name__ == '__main__':
    import time
    t = time.time()
    test_big_solar_with_steps()
    diff = time.time() - t
    print(f'Time: {diff:.2f}')
