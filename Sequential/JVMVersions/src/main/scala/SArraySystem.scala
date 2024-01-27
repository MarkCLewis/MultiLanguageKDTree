class SArraySystem(val numBodies: Int) extends SSystem {
  private val pos = Array.fill(3*numBodies)(0.0)
  private val vel = Array.fill(3*numBodies)(0.0)
  private val radius = Array.fill(numBodies)(0.0)
  private val mass = Array.fill(numBodies)(0.0)

  inline def p(index: Int, dim: Int): Double = pos(3*index + dim)
  inline def v(index: Int, dim: Int): Double = vel(3*index + dim)
  inline def r(index: Int): Double = radius(index)
  inline def m(index: Int): Double = mass(index)

  def init(index: Int, x: Double, y: Double, z: Double, vx: Double, vy: Double, vz: Double, rad: Double, mass: Double): Unit = {
    pos(index*3 + 0) = x
    pos(index*3 + 0) = y
    pos(index*3 + 0) = z
    vel(index*3 + 0) = vx
    vel(index*3 + 0) = vy
    vel(index*3 + 0) = vz
    radius(index) = rad
    this.mass(index) = mass
  }
  
  def incP(index: Int, dim: Int, delta: Double): Unit = pos(index*3 + dim) += delta
  def incV(index: Int, dim: Int, delta: Double): Unit = vel(index*3 + dim) += delta
}
