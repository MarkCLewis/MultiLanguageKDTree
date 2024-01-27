package main

import (
	"fmt"
	"math"
	"testing"
)

func recur_test_tree_struct(
	t *testing.T,
	node int,
	nodes []KDTree,
	particles []Particle,
	min *[3]float64,
	max *[3]float64,
) {
	if nodes[node].num_parts > 0 {
		for index := 0; index < nodes[node].num_parts; index++ {
			i := nodes[node].particles[index]
			for dim := 0; dim < 2; dim++ {
				if particles[i].p[dim] <= min[dim] {
					fmt.Printf("Particle dim %d is below min. i=%d p=%e min=%e\n",
						dim,
						i,
						particles[i].p[dim],
						min[dim])
				}
				if particles[i].p[dim] > max[dim] {
					fmt.Printf("Particle dim %d is above max. i=%d p=%e max=%e\n",
						dim,
						i,
						particles[i].p[dim],
						max[dim])
				}
			}
		}
	} else {
		split_dim := nodes[node].split_dim
		tmin := min[split_dim]
		tmax := max[split_dim]
		max[split_dim] = nodes[node].split_val
		recur_test_tree_struct(t, nodes[node].left, nodes, particles, min, max)
		max[split_dim] = tmax
		min[split_dim] = nodes[node].split_val
		recur_test_tree_struct(t, nodes[node].right, nodes, particles, min, max)
		min[split_dim] = tmin
	}
}

func Test_single_node(t *testing.T) {
	parts := Two_bodies()
	node_vec := allocate_node_vec(len(parts))
	indices := make([]int, len(parts))
	for i := 0; i < len(parts); i++ {
		indices[i] = i
	}
	build_tree(indices, 0, len(parts), parts, 0, &node_vec)
	if node_vec[0].num_parts != len(parts) {
		t.Fatalf("Not all particle in single node. %d %d", node_vec[0].num_parts, len(parts))
	}
}

func Test_two_particle_sim(t *testing.T) {
	parts := Two_bodies()
	indices := make([]int, len(parts))
	for i := 0; i < len(parts); i++ {
		indices[i] = i
	}
	dt := math.Pi / 1000
	Simple_sim(parts, dt, 1000)
	fmt.Printf("1000 steps: %e %e, %e %e\n", parts[0].p[0], parts[0].p[1], parts[1].p[0], parts[1].p[1])
	Simple_sim(parts, dt, 1000)
	fmt.Printf("2000 steps: %e %e, %e %e\n", parts[0].p[0], parts[0].p[1], parts[1].p[0], parts[1].p[1])
}

func Test_two_leaves(t *testing.T) {
	parts := circular_orbits(11)
	node_vec := allocate_node_vec(len(parts))
	indices := make([]int, len(parts))
	for i := 0; i < len(parts); i++ {
		indices[i] = i
	}
	build_tree(indices, 0, len(parts), parts, 0, &node_vec)
	recur_test_tree_struct(
		t,
		0,
		node_vec,
		parts,
		&[3]float64{-1e100, -1e100, -1e100},
		&[3]float64{1e100, 1e100, 1e100},
	)
	if node_vec[0].num_parts != 0 {
		t.Fatalf("Root has %d particle.", node_vec[0].num_parts)
	}
	if node_vec[1].num_parts+node_vec[2].num_parts != len(parts) {
		t.Fatalf("Two leaves don't contain all particles. %d + %d != %d", node_vec[1].num_parts, node_vec[2].num_parts, len(parts))
	}
}

func Test_big_solar(t *testing.T) {
	parts := circular_orbits(50)
	node_vec := allocate_node_vec(len(parts))
	indices := make([]int, len(parts))
	for i := 0; i < len(parts); i++ {
		indices[i] = i
	}
	build_tree(indices, 0, len(parts), parts, 0, &node_vec)
	recur_test_tree_struct(
		t,
		0,
		node_vec,
		parts,
		&[3]float64{-1e100, -1e100, -1e100},
		&[3]float64{1e100, 1e100, 1e100},
	)
}

func Test_big_solar_with_steps(t *testing.T) {
	parts := circular_orbits(50)
	Simple_sim(parts, 1e-3, 10)

	node_vec := allocate_node_vec(len(parts))
	indices := make([]int, len(parts))
	for i := 0; i < len(parts); i++ {
		indices[i] = i
	}
	build_tree(indices, 0, len(parts), parts, 0, &node_vec)
	recur_test_tree_struct(
		t,
		0,
		node_vec,
		parts,
		&[3]float64{-1e100, -1e100, -1e100},
		&[3]float64{1e100, 1e100, 1e100},
	)
}
