#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "kdtree.h"
#include "particle.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

const size_t NEGS[MAX_PARTS] = {0, 0, 0, 0,
                                0, 0, 0}; // TODO - give maximum value here

KDTree_array_t allocate_node_vec(size_t num_parts) {
  size_t num_nodes = 2 * (num_parts / (MAX_PARTS - 1) + 1);
  KDTree_array_t ret = new_kdtree_array_t(num_nodes);
  return ret;
}

// Returns the index of the last Node used in the construction.
size_t build_tree(size_t_array_t *indices, size_t start, size_t end,
                  const Particle_array_t *particles, size_t cur_node,
                  KDTree_array_t *nodes) {
  size_t np = end - start;
  if (np <= MAX_PARTS) {
    if (cur_node >= nodes->size) {
      KDTree_resize(nodes, cur_node + 1);
    }
    nodes->ptr[cur_node].num_parts = np;
    for (size_t i = 0; i < np; ++i) {
      nodes->ptr[cur_node].particles[i] = indices->ptr[start + i];
    }
    return cur_node;
  } else {
    // Pick split dim and value
    double min_arr[3] = {1e100, 1e100, 1e100};
    double max_arr[3] = {-1e100, -1e100, -1e100};
    double m = 0.0;
    double cm[3] = {0.0, 0.0, 0.0};
    for (size_t i = start; i < end; ++i) {
      m += particles->ptr[indices->ptr[i]].m;
      cm[0] += particles->ptr[indices->ptr[i]].m *
               particles->ptr[indices->ptr[i]].p[0];
      cm[1] += particles->ptr[indices->ptr[i]].m *
               particles->ptr[indices->ptr[i]].p[1];
      cm[2] += particles->ptr[indices->ptr[i]].m *
               particles->ptr[indices->ptr[i]].p[2];

      min_arr[0] = MIN(min_arr[0], particles->ptr[indices->ptr[i]].p[0]);
      min_arr[1] = MIN(min_arr[1], particles->ptr[indices->ptr[i]].p[1]);
      min_arr[2] = MIN(min_arr[2], particles->ptr[indices->ptr[i]].p[2]);
      max_arr[0] = MAX(max_arr[0], particles->ptr[indices->ptr[i]].p[0]);
      max_arr[1] = MAX(max_arr[1], particles->ptr[indices->ptr[i]].p[1]);
      max_arr[2] = MAX(max_arr[2], particles->ptr[indices->ptr[i]].p[2]);
    }
    cm[0] /= m;
    cm[1] /= m;
    cm[2] /= m;
    size_t split_dim = 0;
    if (max_arr[1] - min_arr[1] > max_arr[split_dim] - min_arr[split_dim]) {
      split_dim = 1;
    }
    if (max_arr[2] - min_arr[2] > max_arr[split_dim] - min_arr[split_dim]) {
      split_dim = 2;
    }
    double size = max_arr[split_dim] - min_arr[split_dim];

    // Partition particles on split_dim
    size_t mid = (start + end) / 2;
    size_t s = start;
    size_t e = end;
    while (s + 1 < e) {
      double r;
      do {
        r = (double)rand() / RAND_MAX;
      } while (r ==
               1); // there is a tiny chance r = 1, but we want the range
                   // [s,e) (excluding e), and if r = 1 we would produce pivot=e

      size_t pivot = (size_t)(r * (e - s)) + s;

      size_t c = indices->ptr[s];
      indices->ptr[s] = indices->ptr[pivot];
      indices->ptr[pivot] = c;

      size_t low = s + 1;
      size_t high = e - 1;
      while (low <= high) {
        if (particles->ptr[indices->ptr[low]].p[split_dim] <
            particles->ptr[indices->ptr[s]].p[split_dim]) {
          low += 1;
        } else {
          size_t c = indices->ptr[low];
          indices->ptr[low] = indices->ptr[high];
          indices->ptr[high] = c;
          high -= 1;
        }
      }

      size_t c2 = indices->ptr[s];
      indices->ptr[s] = indices->ptr[high];
      indices->ptr[high] = c2;

      if (high < mid) {
        s = high + 1;
      } else if (high > mid) {
        e = high;
      } else {
        s = e;
      }
    }
    double split_val = particles->ptr[indices->ptr[mid]].p[split_dim];

    // Recurse on children and build this node.
    size_t left =
        build_tree(indices, start, mid, particles, cur_node + 1, nodes);
    size_t right = build_tree(indices, mid, end, particles, left + 1, nodes);

    if (cur_node >= nodes->size) {
      KDTree_resize(nodes, cur_node + 1);
    }
    nodes->ptr[cur_node].num_parts = 0;
    nodes->ptr[cur_node].split_dim = split_dim;
    nodes->ptr[cur_node].split_val = split_val;
    nodes->ptr[cur_node].m = m;
    nodes->ptr[cur_node].cm[0] = cm[0];
    nodes->ptr[cur_node].cm[1] = cm[1];
    nodes->ptr[cur_node].cm[2] = cm[2];
    nodes->ptr[cur_node].size = size;
    nodes->ptr[cur_node].left = cur_node + 1;
    nodes->ptr[cur_node].right = left + 1;

    return right;
  }
}

void calc_pp_accel(const Particle *pi, const Particle *pj, double acc[3]) {
  double dp[3] = {pi->p[0] - pj->p[0], pi->p[1] - pj->p[1],
                  pi->p[2] - pj->p[2]};
  double dist = sqrt(dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2]);
  double magi = -pj->m / (dist * dist * dist);
  acc[0] += dp[0] * magi;
  acc[1] += dp[1] * magi;
  acc[2] += dp[2] * magi;
}

void accel_recur(size_t cur_node, size_t p, const Particle_array_t *particles,
                 const KDTree_array_t *nodes, double acc[3]) {
  // println!("accel {}", cur_node);
  if (nodes->ptr[cur_node].num_parts > 0) {
    for (size_t i = 0; i < nodes->ptr[cur_node].num_parts; ++i) {
      if (nodes->ptr[cur_node].particles[i] != p) {
        calc_pp_accel(particles->ptr + p,
                      particles->ptr + nodes->ptr[cur_node].particles[i], acc);
      }
    }
  } else {
    double dp[3];
    dp[0] = particles->ptr[p].p[0] - nodes->ptr[cur_node].cm[0];
    dp[1] = particles->ptr[p].p[1] - nodes->ptr[cur_node].cm[1];
    dp[2] = particles->ptr[p].p[2] - nodes->ptr[cur_node].cm[2];
    double dist_sqr = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];
    // println!("dist = {}, size = {}", dist, nodes[cur_node].size);
    if (nodes->ptr[cur_node].size * nodes->ptr[cur_node].size <
        THETA * THETA * dist_sqr) {
      double dist = sqrt(dist_sqr);
      double magi = -nodes->ptr[cur_node].m / (dist_sqr * dist);
      acc[0] += dp[0] * magi;
      acc[1] += dp[1] * magi;
      acc[2] += dp[2] * magi;
    } else {
      accel_recur(nodes->ptr[cur_node].left, p, particles, nodes, acc);
      accel_recur(nodes->ptr[cur_node].right, p, particles, nodes, acc);
    }
  }
}

void calc_accel(size_t p, const Particle_array_t *particles,
                const KDTree_array_t *nodes, double acc[3]) {
  accel_recur(0, p, particles, nodes, acc);
}

void print_tree(int step, const KDTree_array_t *tree,
                const Particle_array_t *particles) {

  if (step > 9999999) {
    fprintf(stderr, "Step too big!\n");
    exit(EXIT_FAILURE);
  }
  char name[8 + 7];
  snprintf(name, 8 + 7, "tree%d.txt", step);

  FILE *file = fopen(name, "w");
  if (file == NULL) {
    fprintf(stderr, "File couldn't be opened!\n");
    exit(EXIT_FAILURE);
  }

  fprintf(file, "%lu\n", particles->size);

  for (size_t i_iter = 0; i_iter < tree->size; ++i_iter) {
    KDTree *n = tree->ptr + i_iter;
    if (n->num_parts > 0) {
      fprintf(file, "L %lu\n", n->num_parts);

      for (size_t i = 0; i < n->num_parts; ++i) {
        size_t p = n->particles[i];
        fprintf(file, "%f %f %f\n", particles->ptr[p].p[0],
                particles->ptr[p].p[1], particles->ptr[p].p[2]);
      }
    } else {
      fprintf(file, "I %lu %f %lu %lu\n", n->split_dim, n->split_val, n->left,
              n->right);
    }
  }

  fclose(file);
}

void simple_sim(Particle_array_t *bodies, double dt, int steps) {
  vect3_array_t acc = new_vect3_array_t(bodies->size);
  for (size_t i = 0; i < bodies->size; ++i) {
    vect3 a = {0.0, 0.0, 0.0};
    acc.ptr[i] = a;
  }
  KDTree_array_t tree = allocate_node_vec(bodies->size);
  size_t_array_t indices = new_range(0, bodies->size);

  for (int step = 0; step < steps; ++step) {
    for (size_t i = 0; i < bodies->size; ++i) {
      indices.ptr[i] = i;
    }

    build_tree(&indices, 0, bodies->size, bodies, 0, &tree);
    if (step % 10 == 0) {
      print_tree(step, &tree, bodies);
    }
    for (size_t i = 0; i < bodies->size; ++i) {
      calc_accel(i, bodies, &tree, acc.ptr[i].v);
    }
    for (size_t i = 0; i < bodies->size; ++i) {
      bodies->ptr[i].v[0] += dt * acc.ptr[i].v[0];
      bodies->ptr[i].v[1] += dt * acc.ptr[i].v[1];
      bodies->ptr[i].v[2] += dt * acc.ptr[i].v[2];
      bodies->ptr[i].p[0] += dt * bodies->ptr[i].v[0];
      bodies->ptr[i].p[1] += dt * bodies->ptr[i].v[1];
      bodies->ptr[i].p[2] += dt * bodies->ptr[i].v[2];
      acc.ptr[i].v[0] = 0.0;
      acc.ptr[i].v[1] = 0.0;
      acc.ptr[i].v[2] = 0.0;
    }
  }

  FREE_ARRAY(acc);
  FREE_ARRAY(tree);
  FREE_ARRAY(indices);
}

vect3_array_t new_vect3_array_t(size_t elem_count) {
  vect3_array_t a = {elem_count, (vect3 *)calloc(elem_count, sizeof(vect3))};

  if (a.ptr == NULL) {
    fprintf(stderr, "calloc failed: it returned NULL\n");
    exit(EXIT_FAILURE);
  }

  return a;
}

KDTree_array_t new_kdtree_array_t(size_t elem_count) {
  KDTree_array_t a = {elem_count, (KDTree *)calloc(elem_count, sizeof(KDTree))};

  if (a.ptr == NULL) {
    fprintf(stderr, "calloc failed: it returned NULL\n");
    exit(EXIT_FAILURE);
  }

  return a;
}

void KDTree_resize(KDTree_array_t *arr, size_t new_elem_count) {
  if (new_elem_count == arr->size) {
    // do nothing
  } else if (new_elem_count > arr->size) {
    arr->ptr = (KDTree *)realloc(arr->ptr, sizeof(KDTree) * new_elem_count);
    arr->size = new_elem_count;
  } else {
    // new_elem_count < arr.size
    arr->ptr = (KDTree *)realloc(arr->ptr, sizeof(KDTree) * new_elem_count);
    arr->size = new_elem_count;
  }
}
