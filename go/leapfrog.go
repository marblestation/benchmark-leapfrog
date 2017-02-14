package main

import "fmt"
import "math"

// this is a comment

func main() {
    const n_particles int = 2 
    var time float64 = 0
    var time_step float64 = 0.08
    var half_time_step float64 = time_step/2.
    var time_limit float64 = 365.25 * 1e6

    // Arrays (initialized to zero by default)
    //var x [n_particles][3] float64
    //var v [n_particles][3] float64
    //var a [n_particles][3] float64
    //var m [n_particles] float64

	//// Create slices and not arrays, since arrays are passed by copy to func
	// Initialize a ten length slice of empty slices
	x := make([][]float64, n_particles)
	// Initialize those n empty slices
	for i := 0; i < n_particles; i++ {
		x[i] = make([]float64, 3) // initialized to zero by default
	}
	v := make([][]float64, n_particles)
	for i := 0; i < n_particles; i++ {
		v[i] = make([]float64, 3) // initialized to zero by default
	}
	a := make([][]float64, n_particles)
	for i := 0; i < n_particles; i++ {
		a[i] = make([]float64, 3) // initialized to zero by default
	}
	m := make([]float64, n_particles)

    for time <= time_limit {
        integrator_leapfrog_part1(n_particles, x, v, half_time_step)
        time += half_time_step
        gravity_calculate_acceleration(n_particles, m, x, a)
        integrator_leapfrog_part2(n_particles, x, v, a, time_step, half_time_step)
        time += half_time_step
    }
    fmt.Println("Positions:", x)
}

func integrator_leapfrog_part1(n_particles int, x [][]float64, v [][]float64, half_time_step float64) {
	for i := 0; i<n_particles; i++ {
		x[i][0]  += half_time_step * v[i][0]
		x[i][1]  += half_time_step * v[i][1]
		x[i][2]  += half_time_step * v[i][2]
	}
}

func integrator_leapfrog_part2(n_particles int, x [][]float64, v [][]float64, a [][]float64, time_step float64, half_time_step float64) {
	for i := 0; i<n_particles; i++ {
		v[i][0] += time_step * a[i][0]
		v[i][1] += time_step * a[i][1]
		v[i][2] += time_step * a[i][2]
		x[i][0]  += half_time_step * v[i][0]
		x[i][1]  += half_time_step * v[i][1]
		x[i][2]  += half_time_step * v[i][2]
	}
}

func gravity_calculate_acceleration(n_particles int, m []float64, x [][]float64, a [][]float64) {
    G := 6.6742367e-11; // m^3.kg^-1.s^-2
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

