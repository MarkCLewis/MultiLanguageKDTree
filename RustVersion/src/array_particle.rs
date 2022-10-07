

pub struct Particle {
  pub p: [f64; 3],
  pub v: [f64; 3],
  pub r: f64,
  pub m: f64
}

pub fn two_bodies() -> Vec<Particle> {
  let mut bodies = Vec::new();
  bodies.push(Particle { p: [0.0, 0.0, 0.0], 
                         v: [0.0, 0.0, 0.0], r: 1.0, m: 1.0 });
  bodies.push(Particle { p: [1.0, 0.0, 0.0], 
                         v: [0.0, 1.0, 0.0], r: 1e-4, m: 1e-20 });
  bodies
}

pub fn circular_orbits(n: usize) -> Vec<Particle> {
  let mut particle_buf = vec![];
    particle_buf.push(Particle {
      p: [0.0, 0.0, 0.0],
      v: [0.0, 0.0, 0.0],
      r: 0.00465047,
      m: 1.0,
    });

    for i in 0..n {
        let d = 0.1 + ((i as f64) * 5.0 / (n as f64));
        let v = f64::sqrt(1.0 / d);
        let theta = fastrand::f64() * 6.28;
        let x = d * f64::cos(theta);
        let y = d * f64::sin(theta);
        let vx = -v * f64::sin(theta);
        let vy = v * f64::cos(theta);
        particle_buf.push(Particle {
            p: [x, y, 0.0],
            v: [vx, vy, 0.0],
            m: 1e-14,
            r: 1e-7,
        });
    }
    particle_buf
}

pub fn distance_sqr(x1: &[f64; 3], x2: &[f64; 3]) -> f64 {
  let dx = x1[0] - x2[0];
  let dy = x1[1] - x2[1];
  let dz = x1[2] - x2[2];
  dx*dx + dy*dy + dz*dz
}

pub fn distance(x1: &[f64; 3], x2: &[f64; 3]) -> f64 {
  f64::sqrt(distance_sqr(x1, x2))
}

// fn calc_accel(i: usize, j: usize, pi: &Particle, pj: &Particle, acc: &mut Vec<[f64; 3]>) {
//   let dp = pi.p - pj.p;
//   let dp2 = dp * dp;
//   let dist = f64::sqrt(dp2.reduce_sum());
//   let magi = -pj.m / (dist*dist*dist);
//   acc[i] += dp * magi;
//   let magj = pi.m / (dist*dist*dist);
//   acc[j] += dp * magj;
// }

pub fn calc_pp_accel(pi: &Particle, pj: &Particle) -> [f64; 3] {
  let dx = pi.p[0] - pj.p[0];
  let dy = pi.p[1] - pj.p[1];
  let dz = pi.p[2] - pj.p[2];
  let dp2 = dx * dx + dy * dy + dz * dz;
  let dist = f64::sqrt(dp2);
  let magi = -pj.m / (dist*dist*dist);
//   println!("magi={}", magi[0]);
  [magi*dx, magi*dy, magi*dz]
}

pub fn calc_cm_accel(pi: &Particle, m: f64, cm: [f64; 3]) -> [f64; 3] {
  let dx = pi.p[0] - cm[0];
  let dy = pi.p[1] - cm[1];
  let dz = pi.p[2] - cm[2];
  let dp2 = dx * dx + dy * dy + dz * dz;
  let dist = f64::sqrt(dp2);
  let magi =-m / (dist*dp2);
  [dx * magi, dy * magi, dz * magi]
}