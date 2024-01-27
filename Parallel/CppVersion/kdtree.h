#ifndef KDTREE
#define KDTREE

#include<vector>

#include "particle.h"

const size_t MAX_PARTS = 7;
const double THETA = 0.3;
const size_t NEGS[MAX_PARTS] = {}; // TODO - give maximum value here

struct KDTree {
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
};

void simple_sim(vector<Particle> &bodies, double dt, int steps);

#endif