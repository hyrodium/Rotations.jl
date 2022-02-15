function Random.rand(rng::AbstractRNG, ::Random.SamplerType{R}) where R <: Union{<:Rotation{2},<:RotMatrix{2}}
    # Union{<:Rotation{2},<:RotMatrix{2}} seems able to be replaced with Rotation{2},
    # but it avoids the method for general dimensional RotMatrix.
    T = eltype(R)
    if T == Any
        T = Float64
    end

    R(2π * rand(rng, T))
end

# A random rotation can be obtained easily with unit quaternions
# The unit sphere in R⁴ parameterizes quaternion rotations according to the
# Haar measure of SO(3) - see e.g. http://math.stackexchange.com/questions/184086/uniform-distributions-on-the-space-of-rotations-in-3d
function Random.rand(rng::AbstractRNG, ::Random.SamplerType{R}) where R <: Union{<:Rotation{3},<:RotMatrix{3}}
    # Union{<:Rotation{3},<:RotMatrix{3}} seems able to be replaced with Rotation{3},
    # but it avoids the method for general dimensional RotMatrix.
    T = eltype(R)
    if T == Any
        T = Float64
    end

    q = QuatRotation(randn(rng, T), randn(rng, T), randn(rng, T), randn(rng, T))
    return R(q)
end

# A random rotation can be obtained via random matrix and nearest_rotation.
# This is slower than the implementation for 2 or 3 dimensions, but it is possible to create a random matrix on SO(n)
function Random.rand(rng::AbstractRNG, ::Random.SamplerType{R}) where R <: RotMatrix{N} where N
    T = eltype(R)
    if T == Any
        T = Float64
    end

    m = randn(rng, SMatrix{N,N,T})
    return nearest_rotation(m)
end

function Random.rand(rng::AbstractRNG, ::Random.SamplerType{R}) where R <: Union{RotX,RotY,RotZ}
    T = eltype(R)
    if T == Any
        T = Float64
    end

    return R(2π*rand(rng, T))
end

function Random.rand(rng::AbstractRNG, ::Random.SamplerType{R}) where R <: Union{RotXY,RotYZ,RotZX, RotXZ, RotYX, RotZY}
    T = eltype(R)
    if T == Any
        T = Float64
    end

    # Not really sure what this distribution is, but it's also not clear what
    # it should be! rand(RotXY) *is* invariant to pre-rotations by a RotX and
    # post-rotations by a RotY...
    return R(2π*rand(rng, T), 2π*rand(rng, T))
end

function Random.rand(rng::AbstractRNG, ::Random.SamplerType{RP}) where RP <: MRP
    RP(rand(rng, QuatRotation))
end

@inline function Random.rand(rng::AbstractRNG, ::Random.SamplerType{RP}) where RP <: RodriguesParam
    RP(rand(rng, QuatRotation))
end

function Random.rand(rng::AbstractRNG, ::Random.SamplerType{<:QuatRotation{T}}) where T
    _normalize(QuatRotation{T}(randn(rng,T), randn(rng,T), randn(rng,T), randn(rng,T)))
end