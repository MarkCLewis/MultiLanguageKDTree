#ifndef PARTICLE
#define PARTICLE

#include <stddef.h>
#include <sys/_types/_size_t.h>

typedef struct {
  double p[3];
  double v[3];
  double r;
  double m;
} Particle;

typedef struct {
  size_t size;
  Particle *ptr;
} Particle_array_t;

Particle_array_t new_array_t(size_t elem_count);

typedef struct {
  size_t size;
  size_t *ptr;
} size_t_array_t;

size_t_array_t new_size_t_array_t(size_t elem_count);

Particle_array_t two_bodies();
Particle_array_t circular_orbits(size_t n);

#endif
