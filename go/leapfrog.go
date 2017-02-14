package main

import (
	"fmt"
	"math"
)

const nParticles = 2

func main() {
	var time float64
	timeStep := 0.08
	halfTimeStep := timeStep / 2.
	timeLimit := 365.25 * 1e6

	// Arrays (initialized to zero by default)
	x := &[nParticles][3]float64{}
	v := &[nParticles][3]float64{}
	a := &[nParticles][3]float64{}
	m := &[nParticles]float64{}

	m[0] = 0.08                 // M_SUN
	m[1] = 3.0e-6               // M_SUN
	x[1][0] = 0.0162            // AU
	x[1][1] = 6.57192058353e-15 // AU
	x[1][2] = 5.74968548652e-16 // AU
	v[1][0] = -1.48427302304e-14
	v[1][1] = 0.0399408809121
	v[1][2] = 0.00349437429104

	for time <= timeLimit {
		integrator_leapfrog_part1(x, v, halfTimeStep)
		time += halfTimeStep
		gravity_calculate_acceleration(m, x, a)
		integrator_leapfrog_part2(x, v, a, timeStep, halfTimeStep)
		time += halfTimeStep
	}
	fmt.Println("Positions:", x)
}

func integrator_leapfrog_part1(x *[nParticles][3]float64, v *[nParticles][3]float64, halfTimeStep float64) {
	for i := 0; i < nParticles; i++ {
		x[i][0] += halfTimeStep * v[i][0]
		x[i][1] += halfTimeStep * v[i][1]
		x[i][2] += halfTimeStep * v[i][2]
	}
}

func integrator_leapfrog_part2(x *[nParticles][3]float64, v *[nParticles][3]float64, a *[nParticles][3]float64, timeStep float64, halfTimeStep float64) {
	for i := 0; i < nParticles; i++ {
		v[i][0] += timeStep * a[i][0]
		v[i][1] += timeStep * a[i][1]
		v[i][2] += timeStep * a[i][2]
		x[i][0] += halfTimeStep * v[i][0]
		x[i][1] += halfTimeStep * v[i][1]
		x[i][2] += halfTimeStep * v[i][2]
	}
}

func gravity_calculate_acceleration(m *[nParticles]float64, x *[nParticles][3]float64, a *[nParticles][3]float64) {
	G := 6.6742367e-11 // m^3.kg^-1.s^-2
	for i := 0; i < nParticles; i++ {
		a[i][0] = 0
		a[i][1] = 0
		a[i][2] = 0
		for j := 0; j < nParticles; j++ {
			if j == i {
				continue
			}
			dx := x[i][0] - x[j][0]
			dy := x[i][1] - x[j][1]
			dz := x[i][2] - x[j][2]
			r := math.Sqrt(dx*dx + dy*dy + dz*dz)
			prefact := -G / (r * r * r) * m[j]
			a[i][0] += prefact * dx
			a[i][1] += prefact * dy
			a[i][2] += prefact * dz
		}
	}
}
