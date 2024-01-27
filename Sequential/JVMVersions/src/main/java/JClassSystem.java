public class JClassSystem implements JSystem {
  static class Particle {
    public Particle(double x, double y, double z, double vx, double vy, double vz, double rad, double mass) {
      p[0] = x;
      p[1] = y;
      p[2] = z;
      v[0] = vx;
      v[1] = vy;
      v[2] = vz;
      r = rad;
      m = mass;
    }
    public double[] p = new double[3];
    public double[] v = new double[3];
    public double r;
    public double m;
  }

  private Particle[] particles;

  public JClassSystem(int nb) {
    particles = new Particle[nb];
  }

  @Override
  public int numBodies() {
    return particles.length;
  }

  @Override
  public double p(int index, int dim) {
    return particles[index].p[dim];
  }

  @Override
  public double v(int index, int dim) {
    return particles[index].v[dim];
  }

  @Override
  public double r(int index) {
    return particles[index].r;
  }

  @Override
  public double m(int index) {
    return particles[index].m;
  }

  @Override
  public void init(int index, double x, double y, double z, double vx, double vy, double vz, double rad, double mass) {
    particles[index] = new Particle(x, y, z, vx, vy, vz, rad, mass);
  }

  @Override
  public void incP(int index, int dim, double delta) {
    particles[index].p[dim] += delta;
  }

  @Override
  public void incV(int index, int dim, double delta) {
    particles[index].v[dim] += delta;
  }
}
