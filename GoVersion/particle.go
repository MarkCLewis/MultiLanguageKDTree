package main

import "math"

type Particle struct {
	p [3]float64
	v [3]float64
	r float64
	m float64
}

func Two_bodies() []Particle {
	bodies := []Particle{
		{[3]float64{0.0, 0.0, 0.0},
			[3]float64{0.0, 0.0, 0.0}, 1.0, 1.0},
		{[3]float64{1.0, 0.0, 0.0},
			[3]float64{0.0, 1.0, 0.0}, 1e-4, 1e-20},
	}
	return bodies
}

// pub fn distance_sqr(x1: &[f64; 3], x2: &[f64; 3]) -> f64 {
//   let dx = x1[0] - x2[0];
//   let dy = x1[1] - x2[1];
//   let dz = x1[2] - x2[2];
//   dx*dx + dy*dy + dz*dz
// }

// pub fn distance(x1: &[f64; 3], x2: &[f64; 3]) -> f64 {
//   f64::sqrt(distance_sqr(x1, x2))
// }

// fn calc_accel(i: usize, j: usize, pi: &Particle, pj: &Particle, acc: &mut Vec<[f64; 3]>) {
//   let dp = pi.p - pj.p;
//   let dp2 = dp * dp;
//   let dist = f64::sqrt(dp2.reduce_sum());
//   let magi = -pj.m / (dist*dist*dist);
//   acc[i] += dp * magi;
//   let magj = pi.m / (dist*dist*dist);
//   acc[j] += dp * magj;
// }

func Calc_pp_accel(pi *Particle, pj *Particle) [3]float64 {
	dx := pi.p[0] - pj.p[0]
	dy := pi.p[1] - pj.p[1]
	dz := pi.p[2] - pj.p[2]
	dp2 := dx*dx + dy*dy + dz*dz
	dist := math.Sqrt(dp2)
	magi := -pj.m / (dist * dist * dist)
	//   println!("magi={}", magi[0]);
	return [3]float64{magi * dx, magi * dy, magi * dz}
}

// pub fn calc_cm_accel(pi: &Particle, m: f64, cm: [f64; 3]) -> [f64; 3] {
//   let dx = pi.p[0] - cm[0];
//   let dy = pi.p[1] - cm[1];
//   let dz = pi.p[2] - cm[2];
//   let dp2 = dx * dx + dy * dy + dz * dz;
//   let dist = f64::sqrt(dp2);
//   let magi =-m / (dist*dp2);
//   [dx * magi, dy * magi, dz * magi]
// }
