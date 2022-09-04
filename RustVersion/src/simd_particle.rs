
use core_simd::*;

pub struct Particle {
  pub p: f64x4,
  pub v: f64x4,
  pub r: f64,
  pub m: f64
}

pub fn two_bodies() -> Vec<Particle> {
  let mut bodies = Vec::new();
  bodies.push(Particle { p: f64x4::splat(0.0), 
                         v: f64x4::splat(0.0), r: 1.0, m: 1.0 });
  bodies.push(Particle { p: f64x4::from_array([1.0, 0.0, 0.0, 0.0]), 
                         v: f64x4::from_array([0.0, 1.0, 0.0, 0.0]), r: 1e-4, m: 1e-20 });
  bodies
}

pub fn circular_orbits(n: usize) -> Vec<Particle> {
  let mut particle_buf = vec![];
    particle_buf.push(Particle {
      p: f64x4::splat(0.0),
      v: f64x4::splat(0.0),
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
            p: f64x4::from_array([x, y, 0.0, 0.0]),
            v: f64x4::from_array([vx, vy, 0.0, 0.0]),
            m: 1e-14,
            r: 1e-7,
        });
    }
    particle_buf
}

pub fn simple_sim(mut bodies: Vec<Particle>, dt: f64) {
  let dt_vec = f64x4::splat(dt);
  let mut acc = Vec::new();
  for _ in 0..bodies.len() { 
      acc.push(f64x4::splat(0.0))
  };
  for step in 1..1000001 {
      for i in 0..bodies.len()-1 {
          for j in i+1..bodies.len() {
              calc_accel(i, j, &bodies[i], &bodies[j], &mut acc);
          }
      }
      for i in 0..bodies.len() { 
          bodies[i].v += dt_vec * acc[i];
          let dp = dt_vec * bodies[i].v;
          bodies[i].p += dp;
          acc[i] = f64x4::splat(0.0);
      }
      if step % 10000 == 0 {
          println!("{} {} {} {} {}", step, bodies[1].p[0], bodies[1].p[1], bodies[1].v[0], bodies[1].v[1]);
      }
  }
}

pub fn distance_sqr(x1: f64x4, x2: f64x4) -> f64 {
  let dp = x1 - x2;
  let dp2 = dp * dp;
  dp2.reduce_sum()
}

pub fn distance(x1: f64x4, x2: f64x4) -> f64 {
  f64::sqrt(distance_sqr(x1, x2))
}

fn calc_accel(i: usize, j: usize, pi: &Particle, pj: &Particle, acc: &mut Vec<f64x4>) {
  let dp = pi.p - pj.p;
  let dp2 = dp * dp;
  let dist = f64::sqrt(dp2.reduce_sum());
  let magi = f64x4::splat(-pj.m / (dist*dist*dist));
  acc[i] += dp * magi;
  let magj = f64x4::splat(pi.m / (dist*dist*dist));
  acc[j] += dp * magj;
}

pub fn calc_pp_accel(pi: &Particle, pj: &Particle) -> f64x4 {
  let dp = pi.p - pj.p;
  let dp2 = dp * dp;
  let dist = f64::sqrt(dp2.reduce_sum());
  let magi = f64x4::splat(-pj.m / (dist*dist*dist));
//   println!("magi={}", magi[0]);
  dp * magi
}

pub fn calc_cm_accel(pi: &Particle, m: f64, cm: f64x4) -> f64x4 {
  let dp = pi.p - cm;
  let dp2 = dp * dp;
  let dist = f64::sqrt(dp2.reduce_sum());
  let magi = f64x4::splat(-m / (dist*dist*dist));
  dp * magi
}