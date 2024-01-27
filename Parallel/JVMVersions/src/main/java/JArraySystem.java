public class JArraySystem implements JSystem {
  private double[] pos;
  private double[] vel;
  private double[] radius;
  private double[] mass;

  public JArraySystem(int nb) {
    pos = new double[3*nb];
    vel = new double[3*nb];
    radius = new double[nb];
    mass = new double[nb];
  }

  @Override
  public int numBodies() {
    return radius.length;
  }

  @Override
  public double p(int index, int dim) {
    return pos[3*index + dim];
  }

  @Override
  public double v(int index, int dim) {
    return vel[3*index + dim];
  }

  @Override
  public double r(int index) {
    return radius[index];
  }

  @Override
  public double m(int index) {
    return mass[index];
  }

  @Override
  public void init(int index, double x, double y, double z, double vx, double vy, double vz, double rad, double mass) {
    pos[3*index + 0] = x;
    pos[3*index + 1] = y;
    pos[3*index + 2] = z;
    vel[3*index + 0] = vx;
    vel[3*index + 1] = vy;
    vel[3*index + 2] = vz;
    radius[index] = rad;
    this.mass[index] = mass;
  }

  @Override
  public void incP(int index, int dim, double delta) {
    pos[3*index + dim] += delta;
  }

  @Override
  public void incV(int index, int dim, double delta) {
    vel[3*index + dim] += delta;
  }
  
}
