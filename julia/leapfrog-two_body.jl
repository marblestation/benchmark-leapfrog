using StaticArrays

function gravitate!(x::AbstractArray{<:Real, 2}, a::AbstractArray{<:Real, 2},
                   m::AbstractArray{<:Real, 1})
    @assert size(x) == size(a)
    @assert size(m)[1] == size(x)[1]
    @assert size(x)[2] == 3

    G = 6.6742367e-11
    fill!(a, 0)
    for i in 1:length(m)
        for j in 1:length(m)
            if i == j
                continue
            end

            @inbounds dx = x[i,1] - x[j,1]
            @inbounds dy = x[i,2] - x[j,2]
            @inbounds dz = x[i,3] - x[j,3]

            r = hypot(dx, dy, dz)

            @inbounds prefactor = -G/r^3 * m[j]

            @inbounds a[i,1] = muladd(prefactor, dx, a[i,1])
            @inbounds a[i,2] = muladd(prefactor, dy, a[i,2])
            @inbounds a[i,3] = muladd(prefactor, dz, a[i,3])
        end
    end
end

function simulate!(x::AbstractArray{<:Real, 2}, v::AbstractArray{<:Real, 2},
                   a::AbstractArray{<:Real, 2}, m::AbstractArray{<:Real, 1}
                   ; dt::Real=0.08, tmax::Real=3.6525e8)
    @assert size(x) == size(a)
    @assert size(x) == size(v)
    @assert size(m)[1] == size(x)[1]
    @assert size(x)[2] == 3

    hdt = dt / 2
    t = 0
    
    while t < tmax
        @. x += hdt * v
        gravitate!(x, a, m)
        @. v += dt * a
        @. x += hdt * v
        t += dt
    end
end
 
function main(;dt::Real=0.08, tmax::Real=3.6525e8)
    
    # Arrays (initialized to zero by default)
    x = zeros(2,3)
    v = zeros(2,3)
    a = zeros(2,3)

    @inbounds x[2,1] = 0.0162
    @inbounds x[2,2] = 6.57192058353e-15
    @inbounds x[2,3] = 5.74968548652e-16

    @inbounds v[2,1] = -1.48427302304e-14
    @inbounds v[2,2] = 0.0399408809121
    @inbounds v[2,3] = 0.00349437429104

    m = [0.08, 3.0e-6]   # M_Sun

    simulate!(x, v, a, m, dt=dt, tmax=tmax)
    
    println("Positions:", x)
end

function static_main(;dt::Real=0.08, tmax::Real=3.6525e8)
    
    # Arrays (initialized to zero by default)
    x = @MMatrix zeros(2,3)
    v = @MMatrix zeros(2,3)
    a = @MMatrix zeros(2,3)

    @inbounds x[2,1] = 0.0162
    @inbounds x[2,2] = 6.57192058353e-15
    @inbounds x[2,3] = 5.74968548652e-16

    @inbounds v[2,1] = -1.48427302304e-14
    @inbounds v[2,2] = 0.0399408809121
    @inbounds v[2,3] = 0.00349437429104

    m = [0.08, 3.0e-6]   # M_Sun
    
    simulate!(x, v, a, m, dt=dt, tmax=tmax)
    
    println("Positions:", x)
end

