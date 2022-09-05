#include<math.h>
#include<random>

#include "particle.h"

vector<Particle> two_bodies() {
  vector<Particle> bodies;
  bodies.push_back(Particle { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 1.0 });
  bodies.push_back(Particle { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 1e-4, 1e-20 });
  return bodies;
}

vector<Particle> circular_orbits(size_t n) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 6.28);
  vector<Particle> particle_buf;
  particle_buf.push_back(Particle { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0.00465047, 1.0 });

  for (size_t i = 0; i < n; ++i) {
      double d = 0.1 + (i * 5.0 / n);
      double v = sqrt(1.0 / d);
      double theta = dis(gen);
      double x = d * cos(theta);
      double y = d * sin(theta);
      double vx = -v * sin(theta);
      double vy = v * cos(theta);
      particle_buf.push_back(Particle {
          { x, y, 0.0 },
          { vx, vy, 0.0 },
          1e-14,
          1e-7
      });
  }
  return particle_buf;
}