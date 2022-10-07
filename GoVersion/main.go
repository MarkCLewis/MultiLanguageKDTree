package main

import (
	"fmt"
	"math"
	"math/rand"
	"os"
	"strconv"
)

func main() {
	n, err := strconv.Atoi(os.Args[1])
	if err != nil {
		fmt.Println("You need to provide a number of particles and number of steps.")
		return
	}
	steps, err := strconv.Atoi(os.Args[2])
	if err != nil {
		fmt.Println("You need to provide a number of steps.")
		return
	}

	dt := 1e-3 // * 2.0 * std::f64::consts::PI;

	Simple_sim(circular_orbits(n), dt, steps)
}

func circular_orbits(n int) []Particle {
	particle_buf := make([]Particle, n)
	particle_buf[0] = Particle{
		[3]float64{0.0, 0.0, 0.0},
		[3]float64{0.0, 0.0, 0.0},
		0.00465047,
		1.0,
	}

	for i := 0; i < n; i++ {
		d := float64(0.1) + (float64(i) * 5.0 / float64(n))
		v := math.Sqrt(1.0 / d)
		theta := rand.Float64() * 6.28
		x := d * math.Cos(theta)
		y := d * math.Sin(theta)
		vx := -v * math.Sin(theta)
		vy := v * math.Cos(theta)
		particle_buf = append(particle_buf, Particle{
			[3]float64{x, y, 0.0},
			[3]float64{vx, vy, 0.0},
			1e-14,
			1e-7,
		})
	}
	return particle_buf
}
