#ifndef PARTICLE
#define PARTICLE

#include<vector>

using std::vector;

struct Particle {
  double p[3];
  double v[3];
  double r;
  double m;
};

vector<Particle> two_bodies();
vector<Particle> circular_orbits(size_t n);

#endif