import numpy as np

from particle import two_bodies, circular_orbits
from kd_tree import System, recur_test_tree_struct, simple_sim, print_tree


def test_single_node():
    p, v, r, m = two_bodies()
    sys = System.from_amount(len(p))
    assert len(sys.nodes) == 2
    sys.build_tree(0, len(p), (p, v, r, m), 0)
    assert sys.nodes[0].num_parts == len(p)


def test_two_leaves():
    p, v, r, m = circular_orbits(11)
    sys = System.from_amount(len(p))
    assert len(sys.nodes) == 4  # In the rust version, this fails also
    sys.build_tree(0, len(p), (p, v, r, m), 0)
    recur_test_tree_struct(
        0,
        sys.nodes,
        (p, v, r, m),
        np.ones(3, dtype=np.float64) * -1e100,
        np.ones(3, dtype=np.float64) * 1e100,
    )
    assert sys.nodes[0].num_parts == 0
    assert sys.nodes[0].split_dim == 0
    assert sys.nodes[1].num_parts == 5
    assert sys.nodes[2].num_parts == 6


def test_big_solar():
    p, v, r, m = circular_orbits(5000)
    sys = System.from_amount(len(p))
    sys.build_tree(0, len(p), (p, v, r, m), 0)
    recur_test_tree_struct(
        0,
        sys.nodes,
        (p, v, r, m),
        np.ones(3, dtype=np.float64) * (-1e100),
        np.ones(3, dtype=np.float64) * (1e100),
    )

    # # TODO: remove again
    # print_tree(0, sys.nodes, parts)


def test_big_solar_with_steps():
    p, v, r, m = circular_orbits(5000)
    simple_sim((p, v, r, m), 1e-3, 10, print_steps=True)

    sys = System.from_amount(len(p))

    sys.build_tree(0, len(p), (p, v, r, m), 0)
    recur_test_tree_struct(
        0,
        sys.nodes,
        (p, v, r, m),
        np.ones(3, dtype=np.float64) * (-1e100),
        np.ones(3, dtype=np.float64) * (1e100),
    )


if __name__ == '__main__':
    import time
    t = time.time()
    test_big_solar_with_steps()
    diff = time.time() - t
    print(f'Time: {diff:.2f}')
