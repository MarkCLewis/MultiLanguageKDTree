#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kdtree.h"
#include "particle.h"

// void assert_eq_int(int actual, int expected) {
//   if (actual != expected) {
//     fprintf(stderr, "Assertion fail, got %d but expected %d\n", actual,
//             expected);
//     abort();
//   }
// }

#define assert(condition, str, ...)                                            \
  if (!(condition)) {                                                          \
    fprintf(stderr, (str), __VA_ARGS__);                                       \
    abort();                                                                   \
  }

void recur_test_tree_struct(size_t node, KDTree_array_t *nodes,
                            Particle_array_t *particles, double min[3],
                            double max[3]) {
  if (nodes->ptr[node].num_parts > 0) {
    for (size_t index = 0; index < nodes->ptr[node].num_parts; ++index) {
      size_t i = nodes->ptr[node].particles[index];
      for (size_t dim = 0; dim < 2; ++dim) {
        assert(particles->ptr[i].p[dim] >= min[dim],
               "Particle dim %lu is below min. i=%lu p=%f min=%f", dim, i,
               particles->ptr[i].p[dim], min[dim]);
        assert(particles->ptr[i].p[dim] < max[dim],
               "Particle dim %lu is above max. i=%lu p=%f max=%f", dim, i,
               particles->ptr[i].p[dim], max[dim]);
      }
    }
  }

  else {
    size_t split_dim = nodes->ptr[node].split_dim;
    double tmin = min[split_dim];
    double tmax = max[split_dim];
    max[split_dim] = nodes->ptr[node].split_val;
    recur_test_tree_struct(nodes->ptr[node].left, nodes, particles, min, max);
    max[split_dim] = tmax;
    min[split_dim] = nodes->ptr[node].split_val;
    recur_test_tree_struct(nodes->ptr[node].right, nodes, particles, min, max);
    min[split_dim] = tmin;
  }
}

void assert_eq_size_t(size_t actual, size_t expected) {
  if (actual != expected) {
    fprintf(stderr, "Assertion fail, got %lu but expected %lu\n", actual,
            expected);
    abort();
  }
}

void single_node() {
  Particle_array_t parts = two_bodies();
  KDTree_array_t node_vec = allocate_node_vec(parts.size);

  assert_eq_size_t(node_vec.size, 2);
  size_t_array_t indices = new_range(0, parts.size);

  build_tree(&indices, 0, parts.size, &parts, 0, &node_vec);
  assert_eq_size_t(node_vec.ptr[0].num_parts, parts.size);

  FREE_ARRAY(parts);
  FREE_ARRAY(node_vec);
  FREE_ARRAY(indices);
}

void two_leaves() {
  Particle_array_t parts = circular_orbits(11);
  KDTree_array_t node_vec = allocate_node_vec(parts.size);
  assert_eq_size_t(node_vec.size, 6);
  size_t_array_t indices = new_range(0, parts.size);

  build_tree(&indices, 0, parts.size, &parts, 0, &node_vec);

  double min[3] = {-1e100, -1e100, -1e100};
  double max[3] = {1e100, 1e100, 1e100};
  recur_test_tree_struct(0, &node_vec, &parts, min, max);
  assert_eq_size_t(node_vec.ptr[0].num_parts, 0);
  assert_eq_size_t(node_vec.ptr[1].num_parts + node_vec.ptr[2].num_parts, 12);

  FREE_ARRAY(parts);
  FREE_ARRAY(node_vec);
  FREE_ARRAY(indices);
}

void big_solar() {
  Particle_array_t parts = circular_orbits(5000);
  KDTree_array_t node_vec = allocate_node_vec(parts.size);
  size_t_array_t indices = new_range(0, parts.size);
  build_tree(&indices, 0, parts.size, &parts, 0, &node_vec);

  double min[3] = {-1e100, -1e100, -1e100};
  double max[3] = {1e100, 1e100, 1e100};
  recur_test_tree_struct(0, &node_vec, &parts, min, max);

  FREE_ARRAY(parts);
  FREE_ARRAY(node_vec);
  FREE_ARRAY(indices);
}

void big_solar_with_steps() {
  Particle_array_t parts = circular_orbits(5000);
  simple_sim(&parts, 1e-3, 10);

  KDTree_array_t node_vec = allocate_node_vec(parts.size);
  size_t_array_t indices = new_range(0, parts.size);
  build_tree(&indices, 0, parts.size, &parts, 0, &node_vec);
  double min[3] = {-1e100, -1e100, -1e100};
  double max[3] = {1e100, 1e100, 1e100};
  recur_test_tree_struct(0, &node_vec, &parts, min, max);

  FREE_ARRAY(parts);
  FREE_ARRAY(node_vec);
  FREE_ARRAY(indices);
}

int main(int argc, char **argv) {
  if (argc != 2 || strcmp(argv[1], "single_node") == 0) {
    fprintf(stderr, "Running test: single_node\n");
    single_node();
  }

  if (argc != 2 || strcmp(argv[1], "two_leaves") == 0) {
    fprintf(stderr, "Running test: two_leaves\n");
    two_leaves();
  }

  if (argc != 2 || strcmp(argv[1], "big_solar") == 0) {
    fprintf(stderr, "Running test: big_solar\n");
    big_solar();
  }

  if (argc != 2 || strcmp(argv[1], "big_solar_with_steps") == 0) {
    fprintf(stderr, "Running test: big_solar_with_steps\n");
    big_solar_with_steps();
  }

  return EXIT_SUCCESS;
}
