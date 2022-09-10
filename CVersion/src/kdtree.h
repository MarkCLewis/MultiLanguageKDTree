#ifndef KDTREE
#define KDTREE

#include "particle.h"
#include <sys/_types/_size_t.h>

#define MAX_PARTS ((size_t)7)
#define THETA ((double)0.3)

typedef struct {
  // For leaves
  size_t num_parts;
  size_t particles[MAX_PARTS];

  // For internal nodes
  size_t split_dim;
  double split_val;
  double m;
  double cm[3];
  double size;
  size_t left;
  size_t right;
} KDTree;

typedef struct {
  double v[3];
} vect3;

void simple_sim(Particle_array_t *bodies, double dt, int steps);

typedef struct {
  size_t size;
  KDTree *ptr;
} KDTree_array_t;

KDTree_array_t new_kdtree_array_t(size_t elem_count);

typedef struct {
  size_t size;
  vect3 *ptr;
} vect3_array_t;

vect3_array_t new_vect3_array_t(size_t elem_count);

void KDTree_resize(KDTree_array_t *arr, size_t new_elem_count);

#endif
