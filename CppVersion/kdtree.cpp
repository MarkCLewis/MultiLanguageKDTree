#include<random>
#include<fstream>
#include<sstream>

#include "kdtree.h"

std::random_device rd;
std::mt19937 gen(rd());


vector<KDTree> allocate_node_vec(size_t num_parts) {
    size_t num_nodes = 2 * (num_parts / (MAX_PARTS-1) + 1);
    vector<KDTree> ret;
    ret.resize(num_nodes);
    return ret;
}

struct vect3 {
  double v[3];
};

// Returns the index of the last Node used in the construction.
size_t build_tree(
    vector<size_t> &indices,
    size_t start,
    size_t end,
    const vector<Particle> &particles,
    size_t cur_node,
    vector<KDTree> &nodes
) {
    size_t np = end - start;
    if (np <= MAX_PARTS) {
        if (cur_node >= nodes.size()) {
            nodes.resize(cur_node + 1);
        }
        nodes[cur_node].num_parts = np;
        for (size_t i = 0; i < np; ++i) {
            nodes[cur_node].particles[i] = indices[start + i];
        }
        return cur_node;
    } else {
        // Pick split dim and value
        double min[3] = {1e100, 1e100, 1e100};
        double max[3] = {-1e100, -1e100, -1e100};
        double m = 0.0;
        double cm[3] = {0.0, 0.0, 0.0};
        for (size_t i = start; i < end; ++i) {
            m += particles[indices[i]].m;
            cm[0] += particles[indices[i]].m * particles[indices[i]].p[0];
            cm[1] += particles[indices[i]].m * particles[indices[i]].p[1];
            cm[2] += particles[indices[i]].m * particles[indices[i]].p[2];
            min[0] = std::min(min[0], particles[indices[i]].p[0]);
            min[1] = std::min(min[1], particles[indices[i]].p[1]);
            min[2] = std::min(min[2], particles[indices[i]].p[2]);
            max[0] = std::max(max[0], particles[indices[i]].p[0]);
            max[1] = std::max(max[1], particles[indices[i]].p[1]);
            max[2] = std::max(max[2], particles[indices[i]].p[2]);
        }
        cm[0] /= m;
        cm[1] /= m;
        cm[2] /= m;
        size_t split_dim = 0;
        if (max[1] - min[1] > max[split_dim] - min[split_dim]) {
            split_dim = 1;
        }
        if (max[2] - min[2] > max[split_dim] - min[split_dim]) {
            split_dim = 2;
        }
        double size = max[split_dim] - min[split_dim];

        // Partition particles on split_dim
        size_t mid = (start + end) / 2;
        size_t s = start;
        size_t e = end;
        while (s + 1 < e) {
            std::uniform_int_distribution<> distr(s, e-1);
            size_t pivot = distr(gen);
            std::swap(indices[s], indices[pivot]);
            auto low = s + 1;
            auto high = e - 1;
            while (low <= high) {
                if (particles[indices[low]].p[split_dim] < particles[indices[s]].p[split_dim]) {
                    low += 1;
                } else {
                    std::swap(indices[low], indices[high]);
                    high -= 1;
                }
            }
            std::swap(indices[s], indices[high]);
            if (high < mid) {
                s = high + 1;
            } else if (high > mid) {
                e = high;
            } else {
                s = e;
            }
        }
        auto split_val = particles[indices[mid]].p[split_dim];

        // Recurse on children and build this node.
        auto left = build_tree(indices, start, mid, particles, cur_node + 1, nodes);
        auto right = build_tree(indices, mid, end, particles, left + 1, nodes);

        if (cur_node >= nodes.size()) {
            nodes.resize(cur_node + 1);
        }
        nodes[cur_node].num_parts = 0;
        nodes[cur_node].split_dim = split_dim;
        nodes[cur_node].split_val = split_val;
        nodes[cur_node].m = m;
        nodes[cur_node].cm[0] = cm[0];
        nodes[cur_node].cm[1] = cm[1];
        nodes[cur_node].cm[2] = cm[2];
        nodes[cur_node].size = size;
        nodes[cur_node].left = cur_node+1;
        nodes[cur_node].right = left+1;

        return right;
    }
}

void calc_pp_accel(const Particle &pi, const Particle &pj, double acc[3]) {
  double dp[3] = { pi.p[0] - pj.p[0], pi.p[1] - pj.p[1], pi.p[2] - pj.p[2] };
  auto dist = sqrt(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]);
  auto magi = -pj.m / (dist*dist*dist);
  acc[0] += dp[0] * magi;
  acc[1] += dp[1] * magi;
  acc[2] += dp[2] * magi;
}

void accel_recur(size_t cur_node, size_t p, const vector<Particle> &particles, const vector<KDTree> &nodes, double acc[3]) {
    // println!("accel {}", cur_node);
    if (nodes[cur_node].num_parts > 0) {
        for (size_t i = 0; i < nodes[cur_node].num_parts; ++i) {
            if (nodes[cur_node].particles[i] != p) {
                calc_pp_accel(particles[p], particles[nodes[cur_node].particles[i]], acc);
            }
        }
    } else {
        double dp[3];
        dp[0] = particles[p].p[0] - nodes[cur_node].cm[0];
        dp[1] = particles[p].p[1] - nodes[cur_node].cm[1];
        dp[2] = particles[p].p[2] - nodes[cur_node].cm[2];
        auto dist_sqr = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2];
        // println!("dist = {}, size = {}", dist, nodes[cur_node].size);
        if (nodes[cur_node].size * nodes[cur_node].size < THETA * THETA * dist_sqr) {
            auto dist = sqrt(dist_sqr);
            auto magi = -nodes[cur_node].m / (dist_sqr*dist);
            acc[0] += dp[0] * magi;
            acc[1] += dp[1] * magi;
            acc[2] += dp[2] * magi;
        } else {
            accel_recur(nodes[cur_node].left, p, particles, nodes, acc);
            accel_recur(nodes[cur_node].right, p, particles, nodes, acc);
        }
    }
}

void calc_accel(size_t p, const vector<Particle> &particles, const vector<KDTree> &nodes, double acc[3]) {
    accel_recur(0, p, particles, nodes, acc);
}

void print_tree(int step, const vector<KDTree> &tree, const vector<Particle> &particles) {
    std::stringstream sstream;
    sstream << "tree" << step << ".txt";
    std::ofstream file(sstream.str());
    
    file << particles.size() << "\n";
    for (const KDTree &n: tree) {
        if (n.num_parts > 0) {
            file << "L " << n.num_parts << "\n";
            for (size_t i = 0; i < n.num_parts; ++i) {
                auto p = n.particles[i];
                file << particles[p].p[0] << " " << particles[p].p[1] << " " << particles[p].p[2] << "\n";
            }
        } else {
            file << "I " << n.split_dim << " " << n.split_val << " " << n.left << " " << n.right << "\n";
        }
    }
}

void simple_sim(vector<Particle> &bodies, double dt, int steps) {
    vector<vect3> acc;
    for (size_t i = 0; i < bodies.size(); ++i) {
        vect3 a = { 0.0, 0.0, 0.0 };
        acc.push_back(a);
    }
    auto tree = allocate_node_vec(bodies.size());
    vector<size_t> indices;
    indices.resize(bodies.size());
    
    for (int step = 0; step < steps; ++step) {
        for (size_t i = 0; i < bodies.size(); ++i) {
            indices[i] = i;
        }
        build_tree(indices, 0, bodies.size(), bodies, 0, tree);
        // if (step % 10 == 0) {
        //     print_tree(step, tree, bodies);
        // }
        for (size_t i = 0; i < bodies.size(); ++i) {
            calc_accel(i, bodies, tree, acc[i].v);
        }
        for (size_t i = 0; i < bodies.size(); ++i) {
            bodies[i].v[0] += dt * acc[i].v[0];
            bodies[i].v[1] += dt * acc[i].v[1];
            bodies[i].v[2] += dt * acc[i].v[2];
            bodies[i].p[0] += dt * bodies[i].v[0];
            bodies[i].p[1] += dt * bodies[i].v[1];
            bodies[i].p[2] += dt * bodies[i].v[2];
            acc[i].v[0] = 0.0;
            acc[i].v[1] = 0.0;
            acc[i].v[2] = 0.0;
        }
    }
}
