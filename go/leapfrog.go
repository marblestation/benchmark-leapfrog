package main

import (
	"fmt"
	"math"
)

const (
	n_particles = 2
	G           = 6.6742367e-11 // m^3.kg^-1.s^-2
)

func main() {
    var time float64 = 0
    var time_step float64 = 0.08
    var half_time_step float64 = time_step/2.
    var time_limit float64 = 365.25 * 1e6

	//// Create slices and not arrays, since arrays are passed by copy to func
	x := &[n_particles][3]float64{
		{0, 0, 0},
		{0.0162, 6.57192058353e-15, 5.74968548652e-16}, // AU
	}
	v := &[n_particles][3]float64{
		{0, 0, 0},
		{-1.48427302304e-14, 0.0399408809121, 0.00349437429104},
	}
	a := &[n_particles][3]float64{}
	m := &[n_particles]float64{0.08, 3.0e-6} // M_SUN

    for time <= time_limit {
        integrator_leapfrog_part1(x, v, half_time_step)
        time += half_time_step
        gravity_calculate_acceleration(m, x, a)
        integrator_leapfrog_part2(x, v, a, time_step, half_time_step)
        time += half_time_step
    }
    fmt.Println("Positions:", x)
}

func integrator_leapfrog_part1(x, v *[n_particles][3]float64, half_time_step float64) {
	for i := 0; i<n_particles; i++ {
		x[i][0]  += half_time_step * v[i][0]
		x[i][1]  += half_time_step * v[i][1]
		x[i][2]  += half_time_step * v[i][2]
	}
}

func integrator_leapfrog_part2(x, v, a *[n_particles][3]float64, time_step, half_time_step float64) {
	for i := 0; i<n_particles; i++ {
		v[i][0] += time_step * a[i][0]
		v[i][1] += time_step * a[i][1]
		v[i][2] += time_step * a[i][2]
		x[i][0]  += half_time_step * v[i][0]
		x[i][1]  += half_time_step * v[i][1]
		x[i][2]  += half_time_step * v[i][2]
	}
}

func gravity_calculate_acceleration(m *[n_particles]float64, x, a *[n_particles][3]float64) {
	for i := 0; i<n_particles; i++ {
		a[i][0] = 0
		a[i][1] = 0
		a[i][2] = 0
		for j := 0; j<n_particles; j++ {
			if j == i {
				continue
            }
            dx := x[i][0] - x[j][0]
            dy := x[i][1] - x[j][1]
            dz := x[i][2] - x[j][2]
            r := math.Sqrt(dx*dx + dy*dy + dz*dz)
            prefact := -G/(r*r*r) * m[j]
            a[i][0] += prefact * dx
            a[i][1] += prefact * dy
            a[i][2] += prefact * dz
        }
	}
}

