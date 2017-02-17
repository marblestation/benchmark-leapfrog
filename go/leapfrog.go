package main

import (
	"fmt"
	"math"
)

const (
	nParticles = 2
	G          = 6.6742367e-11 // m^3.kg^-1.s^-2
)

func main() {
	var t float64
	dt := 0.08
	hdt := dt / 2.
	tmax := 365.25 * 1e6

	// Arrays (initialized to zero by default)
	x := &[nParticles][3]float64{
		{0, 0, 0},
		{0.0162, 6.57192058353e-15, 5.74968548652e-16}, // AU
	}
	v := &[nParticles][3]float64{
		{0, 0, 0},
		{-1.48427302304e-14, 0.0399408809121, 0.00349437429104},
	}
	a := &[nParticles][3]float64{}
	m := &[nParticles]float64{0.08, 3.0e-6} // M_SUN

	for t <= tmax {
		integrator_leapfrog_part1(x, v, hdt)
		gravity_calculate_acceleration(m, x, a)
		integrator_leapfrog_part2(x, v, a, dt, hdt)
		t += dt
	}
	fmt.Println("Positions:", x)
}

func integrator_leapfrog_part1(x, v *[nParticles][3]float64, hdt float64) {
	for i := 0; i < nParticles; i++ {
		x[i][0] += hdt * v[i][0]
		x[i][1] += hdt * v[i][1]
		x[i][2] += hdt * v[i][2]
	}
}

func integrator_leapfrog_part2(x, v, a *[nParticles][3]float64, dt, hdt float64) {
	for i := 0; i < nParticles; i++ {
		v[i][0] += dt * a[i][0]
		v[i][1] += dt * a[i][1]
		v[i][2] += dt * a[i][2]
		x[i][0] += hdt * v[i][0]
		x[i][1] += hdt * v[i][1]
		x[i][2] += hdt * v[i][2]
	}
}

func gravity_calculate_acceleration(m *[nParticles]float64, x, a *[nParticles][3]float64) {
	dx := x[0][0] - x[1][0]
	dy := x[0][1] - x[1][1]
	dz := x[0][2] - x[1][2]
	r := math.Sqrt(dx*dx + dy*dy + dz*dz)
	prefact := -G / (r * r * r)

	pm0 := prefact * m[0]
	pm1 := prefact * m[1]

	a[0][0] = pm1 * dx
	a[0][1] = pm1 * dy
	a[0][2] = pm1 * dz

	a[1][0] = pm0 * -dx
	a[1][1] = pm0 * -dy
	a[1][2] = pm0 * -dz
}
