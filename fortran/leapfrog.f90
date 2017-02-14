
program leapfrog
    implicit none

    integer :: n_particles
    integer :: error
    real, dimension(:), allocatable :: m
    real, dimension(:,:), allocatable :: x, v, a
    real :: time, time_step, time_limit, half_time_step, time_write

    n_particles = 2
    time = 0
    time_write = 0
    time_step = 0.08d0 ! in days
    time_limit = 365.25d0 * 1.d6
    allocate(m(n_particles), stat=error)
    allocate(x(3,n_particles), stat=error)
    allocate(v(3,n_particles), stat=error)
    allocate(a(3,n_particles), stat=error)
    m(1:n_particles) = 0.0d0
    x(1,1:n_particles) = 0.0d0
    x(2,1:n_particles) = 0.0d0
    x(3,1:n_particles) = 0.0d0
    v(1,1:n_particles) = 0.0d0
    v(2,1:n_particles) = 0.0d0
    v(3,1:n_particles) = 0.0d0
    a(1,1:n_particles) = 0.0d0
    a(2,1:n_particles) = 0.0d0
    a(3,1:n_particles) = 0.0d0
    m(1) = 0.08 ! M_SUN
    m(2) = 3.0e-6 ! M_SUN
    x(1, 2) = 0.0162 ! AU
    x(2, 2) = 6.57192058353e-15 ! AU
    x(3, 2) = 5.74968548652e-16 ! AU
    v(1, 2) = -1.48427302304e-14
    v(2, 2) = 0.0399408809121
    v(3, 2) = 0.00349437429104

    half_time_step = 0.5d0*time_step
    do while (time.le.time_limit)
        call integrator_leapfrog_part1(n_particles,x,v,half_time_step)
        time = time + half_time_step
        call gravity_calculate_acceleration(n_particles,m,x,a)
        call integrator_leapfrog_part2(n_particles,x,v,a,time_step,half_time_step)
        time = time + half_time_step
        !if (time.ge.time_write) then
            !write(*,*) time
            !time_write = time_write + 100.*365.25d0
        !endif
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
