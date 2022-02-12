function Random.rand(rng::AbstractRNG, ::Random.SamplerType{R}) where R <: Union{<:Rotation{2},<:RotMatrix{2}}
    # Union{<:Rotation{2},<:RotMatrix{2}} seems meaningless,
    # but it avoids the method for general dimensional RotMatrix
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
    # Union{<:Rotation{3},<:RotMatrix{3}} seems meaningless,
    # but it avoids the method for general dimensional RotMatrix
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

    m = @SMatrix randn(T,N,N)
    return nearest_rotation(m)
end
