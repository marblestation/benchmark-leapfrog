const N_PARTICLES: usize = 2;

fn main() {
    let mut time: f64 = 0.;
    let time_step: f64 = 0.08;
    let half_time_step: f64 = time_step/2.;
    let time_limit: f64 = 365.25 * 1e6;
    let mut x: [[f64; 3]; N_PARTICLES] = [[0.; 3]; N_PARTICLES];
    let mut v: [[f64; 3]; N_PARTICLES] = [[0.; 3]; N_PARTICLES];
    let mut a: [[f64; 3]; N_PARTICLES] = [[0.; 3]; N_PARTICLES];
    let mut m: [f64; N_PARTICLES] = [0.; N_PARTICLES];
    m[0] = 0.08; // M_SUN
    m[1] = 3.0e-6; // M_SUN
    x[1][0] = 0.0162; // AU
    x[1][1] = 6.57192058353e-15; // AU
    x[1][2] = 5.74968548652e-16; // AU
    v[1][0] = -1.48427302304e-14;
    v[1][1] = 0.0399408809121;
    v[1][2] = 0.00349437429104;

    while time <= time_limit {
        integrator_leapfrog_part1(N_PARTICLES, &mut x, &v, half_time_step);
        time += half_time_step;
        gravity_calculate_acceleration(N_PARTICLES, &m, &x, &mut a);
        integrator_leapfrog_part2(N_PARTICLES, &mut x, &mut v, &a, time_step, half_time_step);
        time += half_time_step;
    }
    //println!("Hello, world!");
    println!("{:?}", x)
}

fn integrator_leapfrog_part1(n_particles: usize, x: &mut [[f64; 3]; N_PARTICLES], v: &[[f64; 3]; N_PARTICLES], half_time_step: f64) {
    for i in 0..n_particles {
		x[i][0]  += half_time_step * v[i][0];
		x[i][1]  += half_time_step * v[i][1];
		x[i][2]  += half_time_step * v[i][2];
    }
}

fn integrator_leapfrog_part2(n_particles: usize, x: &mut [[f64; 3]; N_PARTICLES], v: &mut [[f64; 3]; N_PARTICLES], a: &[[f64; 3]; N_PARTICLES], time_step: f64, half_time_step: f64) {
    for i in 0..n_particles {
		v[i][0] += time_step * a[i][0];
		v[i][1] += time_step * a[i][1];
		v[i][2] += time_step * a[i][2];
		x[i][0] += half_time_step * v[i][0];
		x[i][1] += half_time_step * v[i][1];
		x[i][2] += half_time_step * v[i][2];
    }
}

fn gravity_calculate_acceleration(n_particles: usize, m: &[f64; N_PARTICLES], x: &[[f64; 3]; N_PARTICLES], a: &mut [[f64; 3]; N_PARTICLES]) {
    let g = 6.6742367e-11; // m^3.kg^-1.s^-2
    for i in 0..n_particles {
		a[i][0] = 0.;
		a[i][1] = 0.;
		a[i][2] = 0.;
        for j in 0..n_particles {
			if j == i {
				continue;
            }
            let dx = x[i][0] - x[j][0];
            let dy = x[i][1] - x[j][1];
            let dz = x[i][2] - x[j][2];
            let r = (dx*dx + dy*dy + dz*dz).sqrt();
            let prefact = -g/(r*r*r) * m[j];
            a[i][0] += prefact * dx;
            a[i][1] += prefact * dy;
            a[i][2] += prefact * dz;
        }
	}
}
