trait SSystem {
  def numBodies: Int

  def p(index: Int, dim: Int): Double
  def v(index: Int, dim: Int): Double
  def r(index: Int): Double
  def m(index: Int): Double

  def init(index: Int, x: Double, y: Double, z: Double, vx: Double, vy: Double, vz: Double, rad: Double, mass: Double): Unit
  
  def incP(index: Int, dim: Int, delta: Double): Unit
  def incV(index: Int, dim: Int, delta: Double): Unit
}
