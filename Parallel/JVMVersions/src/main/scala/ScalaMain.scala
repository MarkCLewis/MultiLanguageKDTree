
@main def scalaMain(n: Int, steps: Int): Unit = {
  println("Running sim.")

  val dt = 1e-3

  val system = circular_orbits(n)
  SKDTree.simple_sim(system, dt, steps)
}

def circular_orbits(n: Int): SSystem = {
    var system = new SArraySystem(n+1)
    system.init(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00465047, 1.0)
  
    for (i <- 0 until n) {
        val d = 0.1 + (i * 5.0 / n)
        val v = math.sqrt(1.0 / d)
        val theta = util.Random.nextDouble() * 6.28
        val x = d * math.cos(theta)
        val y = d * math.sin(theta)
        val vx = -v * math.sin(theta)
        val vy = v * math.cos(theta)
        system.init(i+1, 
            x, y, 0.0,
            vx, vy, 0.0,
            1e-14,
            1e-7
        )
    }
    system;
  }