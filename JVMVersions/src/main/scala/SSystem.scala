trait SSystem {
  def p(index: Int, dim: Int): Double
  def v(index: Int, dim: Int): Double
  def r(index: Int): Double
  def m(index: Int): Double

  def clearAccelerations(): Unit
  def incAcc(index: Int, dim: Int, delta: Double): Unit


  def setp(index: Int, dim: Int, value: Double): Unit
  def setv(index: Int, dim: Int, value: Double): Unit
  def setr(index: Int, value: Double): Unit
  def setm(index: Int, value: Double): Unit
}
