using StaticArrays

const ndims = 3
const PosVector = SVector{ndims, Float64}
const zerov = zero(SVector{ndims})

@inline function gravitate!(x::Vector{PosVector}, a::Vector{PosVector},
                            m::Vector)
    G = 6.6742367e-11
    fill!(a, zerov)
    @inbounds for i in 1:length(m)
        for j in 1:length(m)
            if i == j
                continue
            end

            d = x[i] - x[j]
            prefactor = -G/norm(d)^3 * m[j]
            a[i] += prefactor*d
        end
    end
end

function main(;dt::Real=0.08, tmax::Real=3.6525e8)
    
    # Arrays (initialized to zero by default)
    x = [zerov,
         @SVector([0.0162, 6.57192058353e-15, 5.74968548652e-16])]

    v = [zerov,
         @SVector([-1.48427302304e-14, 0.0399408809121, 0.00349437429104])]

    a = [zerov, zerov]

    m = [0.08, 3.0e-6]   # M_Sun

    hdt = dt / 2

    @inbounds for k in 1:round(Int, tmax/dt)
        for i in 1:length(m)
            x[i] += hdt*v[i]
        end
        gravitate!(x, a, m)
        for i in 1:length(m)
            v[i] += a[i]*dt
            x[i] += v[i]*hdt
        end
    end
      
    println("Positions:", x)
end

println("Running up to tmax=10 to compile...")
@time main(tmax=1e2)
println("Running full benchmark...")
@time main()

