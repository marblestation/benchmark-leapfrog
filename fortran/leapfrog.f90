program leapfrog
    implicit none

    integer, parameter :: n_particles = 2
    real, dimension(n_particles) :: m
    real, dimension(3,n_particles) :: x, v, a
    real, dimension(3) :: dR
    real :: t, dt, t_end, r2
    real :: time, time_step, time_limit, half_time_step
    real, parameter :: G = 6.6742367e-11 ! m^3.kg^-1.s^-2
    integer :: i

    time = 0
    time_step = 0.08 ! time step, days
    time_limit = 365.25e6 ! days
    ! Set initial conditions
    m(:) = (/0.08,3.0e-6/) ! M_SUN
    x(:,1) = 0.0
    x(:,2) = (/0.0162,6.57192058353e-15,5.74968548652e-16/) ! AU
    v(:,1) = 0.0
    v(:,2) = (/-1.48427302304e-14,0.0399408809121,0.00349437429104/)
    a(:,:) = 0.0

    half_time_step = 0.5d0*time_step
    do while (time.le.time_limit)
        call integrator_leapfrog_part1(n_particles,x,v,half_time_step)
        time = time + half_time_step
        call gravity_calculate_acceleration(n_particles,m,x,a)
        call integrator_leapfrog_part2(n_particles,x,v,a,time_step,half_time_step)
        time = time + half_time_step
    enddo
    write(*,*) x

end program leapfrog

subroutine integrator_leapfrog_part1(n_particles,x,v,half_time_step)
    implicit none

    ! Input/Output
    integer,intent(in) :: n_particles
    real,intent(out) :: x(3,n_particles)
    real,intent(in) :: v(3,n_particles)
    real,intent(in) :: half_time_step

    ! Local
    integer :: i

    do i=1,n_particles
        ! Positions
        x(1,i) = x(1,i) + half_time_step * v(1,i)
        x(2,i) = x(2,i) + half_time_step * v(2,i)
        x(3,i) = x(3,i) + half_time_step * v(3,i)
    enddo
end subroutine integrator_leapfrog_part1

subroutine integrator_leapfrog_part2(n_particles,x,v,a,time_step,half_time_step)
    implicit none

    ! Input/Output
    integer,intent(in) :: n_particles
    real,intent(out) :: x(3,n_particles), v(3,n_particles)
    real,intent(in) :: a(3,n_particles)
    real,intent(in) :: time_step, half_time_step

    ! Local
    integer :: i

    do i=1,n_particles
        ! Velocities
        v(1,i) = v(1,i) + time_step * a(1,i)
        v(2,i) = v(2,i) + time_step * a(2,i)
        v(3,i) = v(3,i) + time_step * a(3,i)

        ! Positions
        x(1,i) = x(1,i) + half_time_step * v(1,i)
        x(2,i) = x(2,i) + half_time_step * v(2,i)
        x(3,i) = x(3,i) + half_time_step * v(3,i)
    enddo
end subroutine integrator_leapfrog_part2

subroutine gravity_calculate_acceleration(n_particles,m,x,a_grav)
    implicit none

    ! Input/Output
    integer,intent(in) :: n_particles
    real,intent(in) :: x(3,n_particles)
    real,intent(in) :: m(n_particles)
    real, intent(out) :: a_grav(3,n_particles)

    ! Local
    integer :: i,j
    real :: dx,dy,dz,rr,prefact,G
    !-------------------------------------------------------------------------

    G = 6.6742367e-11 ! m^3.kg^-1.s^-2
    do i = 1,n_particles
        ! Initialization
        a_grav(1,i) = 0.0d0 
        a_grav(2,i) = 0.0d0 
        a_grav(3,i) = 0.0d0 

        do j = 1,n_particles
            if (i.ne.j) then
                dx = x(1,i) - x(1,j)
                dy = x(2,i) - x(2,j)
                dz = x(3,i) - x(3,j)
                rr = sqrt(dx*dx + dy*dy + dz*dz)
                prefact = -G*m(j)/(rr*rr*rr)

                a_grav(1,i) =  a_grav(1,i) + prefact * dx 
                a_grav(2,i) =  a_grav(2,i) + prefact * dy 
                a_grav(3,i) =  a_grav(3,i) + prefact * dz 
            endif
        enddo
    enddo

    !-------------------------------------------------------------------------
    return
end subroutine gravity_calculate_acceleration
