import scala.collection.mutable;
import Particles.*

class SKDTree(
  // For leaves
  var num_parts: Int = 0,
  val particles: Array[Particle] = Array.fill(SKDTree.MAX_PARTS)(Particle(0)),

  // For internal nodes
  var split_dim: Int = 0,
  var split_val: Double = 0.0,
  var m: Double = 0.0,
  val cm: Array[Double] = Array.fill(3)(0.0),
  var size: Double = 0.0,
  var left: Int = 0,
  var right: Int = 0
)

object SKDTree {
  val MAX_PARTS = 7
  val THETA = 0.3

  def allocate_node_vec(numParts: Int): mutable.Buffer[SKDTree] = {
    val numNodes = 2 * (numParts / (MAX_PARTS - 1) + 1)
    mutable.ArrayBuffer.fill(numNodes)(SKDTree())
  }

  // Returns the index of the last Node used in the construction.
  def build_tree(
      indices: Array[Particle],
      start: Int,
      end: Int,
      system: SSystem,
      cur_node: Int,
      nodes: mutable.Buffer[SKDTree]): Int = {
    given SSystem = system

    val np = end - start
    if (np <= MAX_PARTS) {
      while (cur_node >= nodes.size) {
        nodes += SKDTree()
      }
      nodes(cur_node).num_parts = np
      for (i <- 0 until np) {
        nodes(cur_node).particles(i) = indices(start + i)
      }
      return cur_node
    } else {
      // Pick split dim and value
      val min = Array(1e100, 1e100, 1e100)
      val max = Array(-1e100, -1e100, -1e100)
      var m = 0.0
      val cm = Array(0.0, 0.0, 0.0)
      for (i <- start until end) {
        val pi = indices(i)

        m += pi.m// system.m(indices[i]);
        cm(0) += pi.m * pi.p(0)
        cm(1) += pi.m * pi.p(1)
        cm(2) += pi.m * pi.p(2)
        min(0) = Math.min(min(0), pi.p(0))
        min(1) = Math.min(min(1), pi.p(1))
        min(2) = Math.min(min(2), pi.p(2))
        max(0) = Math.max(max(0), pi.p(0))
        max(1) = Math.max(max(1), pi.p(1))
        max(2) = Math.max(max(2), pi.p(2))
      }
      cm(0) /= m
      cm(1) /= m
      cm(2) /= m
      var split_dim = (0 to 2).maxBy(d => max(d) - min(d))
      val size = max(split_dim) - min(split_dim)

      // Partition particles on split_dim
      val mid = (start + end) / 2
      var s = start
      var e = end
      while (s + 1 < e) {
        val pivot = util.Random.nextInt(e-s) + s
        val swapTmp = indices(s)
        indices(s) = indices(pivot)
        indices(pivot) = swapTmp
        var low = s + 1
        var high = e - 1
        while (low <= high) {
          if (indices(low).p(split_dim) < indices(s).p(split_dim)) {
            low += 1
          } else {
            val swapTmp2 = indices(low)
            indices(low) = indices(high)
            indices(high) = swapTmp2
            high -= 1
          }
        }
        val swapTmp3 = indices(s)
        indices(s) = indices(high)
        indices(high) = swapTmp3
        if (high < mid) {
          s = high + 1
        } else if (high > mid) {
          e = high
        } else {
          s = e
        }
      }
      var split_val = indices(mid).p(split_dim)

      // Recurse on children and build this node.
      var left = build_tree(indices, start, mid, system, cur_node + 1, nodes)
      var right = build_tree(indices, mid, end, system, left + 1, nodes)

      while (cur_node >= nodes.size) {
        nodes += SKDTree()
      }
      nodes(cur_node).num_parts = 0;
      nodes(cur_node).split_dim = split_dim;
      nodes(cur_node).split_val = split_val;
      nodes(cur_node).m = m;
      nodes(cur_node).cm(0) = cm(0);
      nodes(cur_node).cm(1) = cm(1);
      nodes(cur_node).cm(2) = cm(2);
      nodes(cur_node).size = size;
      nodes(cur_node).left = cur_node + 1;
      nodes(cur_node).right = left + 1;

      right;
    }
  }

  def calc_pp_accel(system: SSystem, pi: Particle, pj: Particle, acc: Array[Double]): Unit = {
    given SSystem = system
    val dx = pi.p(0) - pj.p(0)
    val dy = pi.p(1) - pj.p(1)
    val dz = pi.p(2) - pj.p(2)
    val dist = math.sqrt(dx * dx + dy * dy + dz * dz)
    val magi = -pj.m / (dist * dist * dist)
    // if(j == 0 && (1.0 - dist).abs < 1e-2) {
    //   System.out.println("magi = " + magi);
    // }
    acc(0) += dx * magi
    acc(1) += dy * magi
    acc(2) += dz * magi
  }

  def accel_recur(cur_node: Int, p: Particle, system: SSystem, nodes: mutable.Buffer[SKDTree], acc: Array[Double]): Unit = {
    given SSystem = system
    // println!("accel {}", cur_node);
    if (nodes(cur_node).num_parts > 0) {
      for (i <- 0 until nodes(cur_node).num_parts) {
        if (nodes(cur_node).particles(i) != p) {
          calc_pp_accel(system, p, nodes(cur_node).particles(i), acc)
        }
      }
    } else {
      val dx = p.p(0) - nodes(cur_node).cm(0)
      val dy = p.p(1) - nodes(cur_node).cm(1)
      val dz = p.p(2) - nodes(cur_node).cm(2)
      var dist_sqr = dx * dx + dy * dy + dz * dz
      // println!("dist = {}, size = {}", dist, nodes[cur_node].size);
      if (nodes(cur_node).size * nodes(cur_node).size < THETA * THETA * dist_sqr) {
        var dist = math.sqrt(dist_sqr)
        var magi = -nodes(cur_node).m / (dist_sqr * dist)
        acc(0) += dx * magi
        acc(1) += dy * magi
        acc(2) += dz * magi
      } else {
        accel_recur(nodes(cur_node).left, p, system, nodes, acc)
        accel_recur(nodes(cur_node).right, p, system, nodes, acc)
      }
    }
  }

  def calc_accel(p: Particle, system: SSystem, nodes: mutable.Buffer[SKDTree], acc: Array[Double]): Unit = {
    accel_recur(0, p, system, nodes, acc)
  }

  def print_tree(step: Int, tree: mutable.Buffer[SKDTree], system: SSystem): Unit = {
    given SSystem = system
    val fname = s"tree$step.txt"
    try {
      var pw = new java.io.PrintWriter(fname)

      pw.println(tree.size)
      for (n <- tree) {
        if (n.num_parts > 0) {
          pw.println(s"L ${n.num_parts}")
          for (i <- 0 until n.num_parts) {
            val p = n.particles(i)
            pw.println(s"${p.p(0)} ${p.p(1)} ${p.p(2)}")
          }
        } else {
          pw.println(s"I ${n.split_dim} ${n.split_val} ${n.left} ${n.right}")
        }
      }
      pw.close()
    } catch {
      case ex: java.io.IOException =>
        System.out.println("Exception writing to file.\n")
        ex.printStackTrace()
    }
  }

  def simple_sim(system: SSystem, dt: Double, steps: Int): Unit = {
    given SSystem = system
    val acc = Array.fill(system.numBodies, 3)(0.0)

    val tree = allocate_node_vec(system.numBodies);
    val indices = Array.fill(system.numBodies)(Particle(0));

    for (step <- 0 until steps) {
      println(step)
      for (i <- 0 until system.numBodies) {
        indices(i) = Particle(i)
      }
      val bstart = System.nanoTime()
      build_tree(indices, 0, system.numBodies, system, 0, tree)
      println(s"Build took: ${(System.nanoTime() - bstart)*1e-9}")
      // if (step % 10 == 0) {
      //   print_tree(step, tree, system)
      // }
      val astart = System.nanoTime()
      for (i <- 0 until system.numBodies) {
        calc_accel(Particle(i), system, tree, acc(i))
      }
      println(s"Acc took: ${(System.nanoTime() - bstart)*1e-9}")
      for (i <- 0 until system.numBodies) {
        val pi = Particle(i)
        pi.incV(0, dt * acc(i)(0))
        pi.incV(1, dt * acc(i)(1))
        pi.incV(2, dt * acc(i)(2))
        pi.incP(0, dt * pi.v(0))
        pi.incP(1, dt * pi.v(1))
        pi.incP(2, dt * pi.v(2))
        acc(i)(0) = 0.0
        acc(i)(1) = 0.0
        acc(i)(2) = 0.0
      }
    }
  }
}
