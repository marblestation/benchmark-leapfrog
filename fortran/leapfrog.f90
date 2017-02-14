program leapfrog
    implicit none

    integer, parameter :: n_particles = 2
    real, dimension(n_particles) :: m
    real, dimension(3,n_particles) :: x, v, a
    real, dimension(3) :: dR
    real :: t, dt, t_end, r2
    real, parameter :: G = 6.6742367e-11 ! m^3.kg^-1.s^-2
    integer :: i

    t = 0
    dt = 0.08 ! time step, days
    t_end = 365.25e6 ! days
    ! Set initial conditions
    m(:) = (/0.08,3.0e-6/) ! M_SUN
    x(:,1) = 0.0
    x(:,2) = (/0.0162,6.57192058353e-15,5.74968548652e-16/) ! AU
    v(:,1) = 0.0
    v(:,2) = (/-1.48427302304e-14,0.0399408809121,0.00349437429104/)
    a(:,:) = 0.0

    do while (t <= t_end)
        ! First step of leapfrog
        x = x + 0.5*dt*v
        t = t + 0.5*dt

        ! Compute forces
        dR(:) = x(:,1) - x(:,2)
        r2 = sqrt(sum(dR**2))
        a(:,1) = -G*m(2)/(r2*r2*r2) * dR
        a(:,2) = -a(:,1) ! Newton's third law

        ! Second step of leapfrog
        v = v + dt*a
        x = x + 0.5*dt*v
        t = t + 0.5*dt
    enddo
    write(*,*) x

end program leapfrog
