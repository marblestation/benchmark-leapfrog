using StaticArrays
using LinearAlgebra: norm

function gravitate!(x::Vector{T}, a::Vector{T},m::AbstractVector) where T
    G = 6.6742367e-11
    fill!(a, zero(T))
    for i in eachindex(x)
        for j in eachindex(x)
            if i == j
                continue
            end
            d = x[i] - x[j]
            prefactor = -G/norm(d)^3 * m[j]
            a[i] += prefactor*d
        end
    end
end

function main(
    x::Vector{T},v::Vector{T},m::AbstractVector;
    dt::Real=0.08, tmax::Real=3.6525e8
) where T
    a = similar(x)
    hdt = dt / 2
    for k in 1:round(Int,tmax/dt)
        @. x += hdt * v
        gravitate!(x, a, m)
        @. v += a*dt
        @. x += v*hdt
    end
    println("Positions:", x)
end

println("Running full benchmark...")

# Dimension of particles
T = SVector{3,Float64}

# Initial positions and velocities
x = [zero(T),T(0.0162, 6.57192058353e-15, 5.74968548652e-16)]
v = [zero(T),T(-1.48427302304e-14, 0.0399408809121, 0.00349437429104)]

# Masses
m = @SVector [0.08, 3.0e-6]

# Run simulation
main(x,v,m)



