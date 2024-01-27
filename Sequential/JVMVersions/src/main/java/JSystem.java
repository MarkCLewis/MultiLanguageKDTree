public interface JSystem {
  int numBodies();

  double p(int index, int dim);
  double v(int index, int dim);
  double r(int index);
  double m(int index);

  void init(int index, double x, double y, double z, double vx, double vy, double vz, double rad, double mass);
  
  void incP(int index, int dim, double delta);
  void incV(int index, int dim, double delta);
}
