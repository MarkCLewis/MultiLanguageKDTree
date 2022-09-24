import java.util.random.RandomGenerator;

public class JavaMain {
  public static void main(String[] args) {
    if (args.length < 2) {
      System.out.println("Specify a number of particles and number of steps.");
      return;
    }
    System.out.println("Running sim.");
  
    int n = Integer.parseInt(args[0]);
    int steps = Integer.parseInt(args[1]);
  
    double dt = 1e-3; // * 2.0 * std::f64::consts::PI;
  
    // let start = Instant::now();
    var system = circular_orbits(n);
    JKDTree.simple_sim(system, dt, steps);
    // println!("{}", start.elapsed().as_nanos() as f64 / 1e9);
  }

  static JSystem circular_orbits(int n) {
    var rand = RandomGenerator.getDefault();
    var system = new JArraySystem(n+1);
    system.init(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00465047, 1.0 );
  
    for (int i = 0; i < n; ++i) {
        double d = 0.1 + (i * 5.0 / n);
        double v = Math.sqrt(1.0 / d);
        double theta = rand.nextDouble(0.0, 6.28);
        double x = d * Math.cos(theta);
        double y = d * Math.sin(theta);
        double vx = -v * Math.sin(theta);
        double vy = v * Math.cos(theta);
        system.init(i+1, 
            x, y, 0.0,
            vx, vy, 0.0,
            1e-14,
            1e-7
        );
    }
    return system;
  }
}
