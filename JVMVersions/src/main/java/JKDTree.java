import java.io.IOException;
import java.util.ArrayList;
import java.util.random.RandomGenerator;

public class JKDTree {
  static final int MAX_PARTS = 7;
  static final double THETA = 0.3;
  static final RandomGenerator rand = RandomGenerator.getDefault();

  // For leaves
  private int num_parts;
  private int[] particles = new int[JKDTree.MAX_PARTS];

  // For internal nodes
  private int split_dim;
  private double split_val;
  private double m;
  private double[] cm = new double[3];
  private double size;
  private int left;
  private int right;

  static ArrayList<JKDTree> allocate_node_vec(int num_parts) {
    int num_nodes = 2 * (num_parts / (MAX_PARTS - 1) + 1);
    ArrayList<JKDTree> ret = new ArrayList<>(num_nodes);
    return ret;
  }

  // Returns the index of the last Node used in the construction.
  static int build_tree(
      int[] indices,
      int start,
      int end,
      JSystem system,
      int cur_node,
      ArrayList<JKDTree> nodes) {
    int np = end - start;
    if (np <= MAX_PARTS) {
      while (cur_node >= nodes.size()) {
        nodes.add(new JKDTree());
      }
      nodes.get(cur_node).num_parts = np;
      for (int i = 0; i < np; ++i) {
        nodes.get(cur_node).particles[i] = indices[start + i];
      }
      return cur_node;
    } else {
      // Pick split dim and value
      double[] min = { 1e100, 1e100, 1e100 };
      double[] max = { -1e100, -1e100, -1e100 };
      double m = 0.0;
      double[] cm = { 0.0, 0.0, 0.0 };
      for (int i = start; i < end; ++i) {
        m += system.m(indices[i]);
        cm[0] += system.m(indices[i]) * system.p(indices[i], 0);
        cm[1] += system.m(indices[i]) * system.p(indices[i], 1);
        cm[2] += system.m(indices[i]) * system.p(indices[i], 2);
        min[0] = Math.min(min[0], system.p(indices[i], 0));
        min[1] = Math.min(min[1], system.p(indices[i], 1));
        min[2] = Math.min(min[2], system.p(indices[i], 2));
        max[0] = Math.max(max[0], system.p(indices[i], 0));
        max[1] = Math.max(max[1], system.p(indices[i], 1));
        max[2] = Math.max(max[2], system.p(indices[i], 2));
      }
      cm[0] /= m;
      cm[1] /= m;
      cm[2] /= m;
      int split_dim = 0;
      if (max[1] - min[1] > max[split_dim] - min[split_dim]) {
        split_dim = 1;
      }
      if (max[2] - min[2] > max[split_dim] - min[split_dim]) {
        split_dim = 2;
      }
      double size = max[split_dim] - min[split_dim];

      // Partition particles on split_dim
      int mid = (start + end) / 2;
      int s = start;
      int e = end;
      while (s + 1 < e) {
        int pivot = rand.nextInt(s, e);
        int swapTmp = indices[s];
        indices[s] = indices[pivot];
        indices[pivot] = swapTmp;
        var low = s + 1;
        var high = e - 1;
        while (low <= high) {
          if (system.p(indices[low], split_dim) < system.p(indices[s], split_dim)) {
            low += 1;
          } else {
            int swapTmp2 = indices[low];
            indices[low] = indices[high];
            indices[high] = swapTmp2;
            high -= 1;
          }
        }
        int swapTmp3 = indices[s];
        indices[s] = indices[high];
        indices[high] = swapTmp3;
        if (high < mid) {
          s = high + 1;
        } else if (high > mid) {
          e = high;
        } else {
          s = e;
        }
      }
      var split_val = system.p(indices[mid], split_dim);

      // Recurse on children and build this node.
      var left = build_tree(indices, start, mid, system, cur_node + 1, nodes);
      var right = build_tree(indices, mid, end, system, left + 1, nodes);

      while (cur_node >= nodes.size()) {
        nodes.add(new JKDTree());
      }
      nodes.get(cur_node).num_parts = 0;
      nodes.get(cur_node).split_dim = split_dim;
      nodes.get(cur_node).split_val = split_val;
      nodes.get(cur_node).m = m;
      nodes.get(cur_node).cm[0] = cm[0];
      nodes.get(cur_node).cm[1] = cm[1];
      nodes.get(cur_node).cm[2] = cm[2];
      nodes.get(cur_node).size = size;
      nodes.get(cur_node).left = cur_node + 1;
      nodes.get(cur_node).right = left + 1;

      return right;
    }
  }

  static void calc_pp_accel(JSystem system, int i, int j, double[] acc) {
    double dx = system.p(i, 0) - system.p(j, 0);
    double dy = system.p(i, 1) - system.p(j, 1);
    double dz = system.p(i, 2) - system.p(j, 2);
    var dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
    var magi = -system.m(j) / (dist * dist * dist);
    if (j == 0 && Math.abs(1.0 - dist) < 1e-2) {
      System.out.println("magi = " + magi);
    }
    acc[0] += dx * magi;
    acc[1] += dy * magi;
    acc[2] += dz * magi;
  }

  static void accel_recur(int cur_node, int p, JSystem system, ArrayList<JKDTree> nodes, double[] acc) {
    // println!("accel {}", cur_node);
    if (nodes.get(cur_node).num_parts > 0) {
      for (int i = 0; i < nodes.get(cur_node).num_parts; ++i) {
        if (nodes.get(cur_node).particles[i] != p) {
          calc_pp_accel(system, p, nodes.get(cur_node).particles[i], acc);
        }
      }
    } else {
      double dx = system.p(p, 0) - nodes.get(cur_node).cm[0];
      double dy = system.p(p, 1) - nodes.get(cur_node).cm[1];
      double dz = system.p(p, 2) - nodes.get(cur_node).cm[2];
      var dist_sqr = dx * dx + dy * dy + dz * dz;
      // println!("dist = {}, size = {}", dist, nodes[cur_node].size);
      if (nodes.get(cur_node).size * nodes.get(cur_node).size < THETA * THETA * dist_sqr) {
        var dist = Math.sqrt(dist_sqr);
        var magi = -nodes.get(cur_node).m / (dist_sqr * dist);
        acc[0] += dx * magi;
        acc[1] += dy * magi;
        acc[2] += dz * magi;
      } else {
        accel_recur(nodes.get(cur_node).left, p, system, nodes, acc);
        accel_recur(nodes.get(cur_node).right, p, system, nodes, acc);
      }
    }
  }

  static void calc_accel(int p, JSystem system, ArrayList<JKDTree> nodes, double[] acc) {
    accel_recur(0, p, system, nodes, acc);
  }

  static void print_tree(int step, ArrayList<JKDTree> tree, JSystem system) {
    String fname = "tree" + step + ".txt";
    try {
      var pw = new java.io.PrintWriter(fname);

      pw.println(tree.size());
      for (JKDTree n : tree) {
        if (n.num_parts > 0) {
          pw.println("L " + n.num_parts);
          for (int i = 0; i < n.num_parts; ++i) {
            var p = n.particles[i];
            pw.println(system.p(p, 0) + " " + system.p(p, 1) + " " + system.p(p, 2));
          }
        } else {
          pw.println("I " + n.split_dim + " " + n.split_val + " " + n.left + " " + n.right);
        }
      }
      pw.close();
    } catch (IOException ex) {
      System.out.println("Exception writing to file.\n");
      ex.printStackTrace();
    }
  }

  static void simple_sim(JSystem system, double dt, int steps) {
    double[][] acc = new double[system.numBodies()][3];

    var tree = allocate_node_vec(system.numBodies());
    int[] indices = new int[system.numBodies()];

    for (int step = 0; step < steps; ++step) {
      for (int i = 0; i < system.numBodies(); ++i) {
        indices[i] = i;
      }
      var bstart = System.nanoTime();
      build_tree(indices, 0, system.numBodies(), system, 0, tree);
      System.out.println("Build took: " + (System.nanoTime() - bstart)*1e-9);
      // if (step % 10 == 0) {
      //   print_tree(step, tree, system);
      // }
      var astart = System.nanoTime();
      for (int i = 0; i < system.numBodies(); ++i) {
        calc_accel(i, system, tree, acc[i]);
      }
      System.out.println("Build took: " + (System.nanoTime() - astart)*1e-9);
      for (int i = 0; i < system.numBodies(); ++i) {
        system.incV(i, 0, dt * acc[i][0]);
        system.incV(i, 1, dt * acc[i][1]);
        system.incV(i, 2, dt * acc[i][2]);
        system.incP(i, 0, dt * system.v(i, 0));
        system.incP(i, 1, dt * system.v(i, 1));
        system.incP(i, 2, dt * system.v(i, 2));
        acc[i][0] = 0.0;
        acc[i][1] = 0.0;
        acc[i][2] = 0.0;
      }
    }
  }

}
