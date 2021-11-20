"""
    abstract type InfinitesimalRotation{N,T} <: StaticMatrix{N,N,T}

An abstract type representing `N`-dimensional infinitesimal rotations. More abstractly, they represent
skew-symmetric real `N`×`N` matrices.
"""
abstract type InfinitesimalRotation{N,T} <: StaticMatrix{N,N,T} end

Base.@pure StaticArrays.Size(::Type{InfinitesimalRotation{N}}) where {N} = Size(N,N)
Base.@pure StaticArrays.Size(::Type{InfinitesimalRotation{N,T}}) where {N,T} = Size(N,N)
Base.@pure StaticArrays.Size(::Type{R}) where {R<:InfinitesimalRotation} = Size(supertype(R))
Base.adjoint(r::InfinitesimalRotation) = -r
Base.transpose(r::InfinitesimalRotation{N,T}) where {N,T<:Real} = -r

# Generate identity-matrix with SMatrix
# Note that zeros(InfinitesimalRotation3,dims...) is not Array{<:InfinitesimalRotation} but Array{<:StaticMatrix{3,3}}
Base.one(::InfinitesimalRotation{N,T}) where {N,T} = one(SMatrix{N,N,T})
Base.one(::Type{InfinitesimalRotation}) = error("The dimension of rotation is not specified.")
Base.one(::Type{<:InfinitesimalRotation{N}}) where N = one(SMatrix{N,N})
Base.one(::Type{<:InfinitesimalRotation{N,T}}) where {N,T} = one(SMatrix{N,N,T})
Base.ones(::Type{R}) where {R<:InfinitesimalRotation} = ones(R, ()) # avoid StaticArray constructor
# Base.ones(::Type{R}, dims::Base.DimOrInd...) where {R<:InfinitesimalRotation} = ones(typeof(one(R)),dims...)
# Base.ones(::Type{R}, dims::NTuple{N, Integer}) where {R<:InfinitesimalRotation, N} = ones(typeof(one(R)),dims)

# Generate identity rotation matrix
Base.zero(r::InfinitesimalRotation) = zero(typeof(r))

# `convert` goes through the constructors, similar to e.g. `Number`
Base.convert(::Type{R}, rot::InfinitesimalRotation{N}) where {N,R<:InfinitesimalRotation{N}} = R(rot)

# InfinitesimalRotation matrices should be orthoginal/unitary. Only the operations we define,
# like multiplication, will stay as InfinitesimalRotations, otherwise users will get an
# SMatrix{3,3} (e.g. rot1 + rot2 -> SMatrix)
Base.@pure StaticArrays.similar_type(::Union{R,Type{R}}) where {R <: InfinitesimalRotation} = SMatrix{size(R)..., eltype(R), prod(size(R))}
Base.@pure StaticArrays.similar_type(::Union{R,Type{R}}, ::Type{T}) where {R <: InfinitesimalRotation, T} = SMatrix{size(R)..., T, prod(size(R))}

# Useful for converting arrays of rotations to another rotation eltype, for instance.
# Only works because parameters of all the rotations are of a similar form
# Would need to be more sophisticated if we have arbitrary dimensions, etc
@inline function Base.promote_op(::Type{R1}, ::Type{R2}) where {R1 <: InfinitesimalRotation, R2 <: InfinitesimalRotation}
    size(R1) == size(R2) || throw(DimensionMismatch("cannot promote rotations of $(size(R1)[1]) and $(size(R2)[1]) dimensions"))
    if isleaftype(R1)
        return R1
    else
        return R1{eltype(R2)}
    end
end

################################################################################
################################################################################
"""
    struct InfinitesimalRotMatrix{N,T} <: InfinitesimalRotation{N,T}

A statically-sized, N×N skew-symmetric matrix.

Note: the skew-symmetricity of the input matrix is *not* checked by the constructor.
"""
struct InfinitesimalRotMatrix{N,T,L} <: InfinitesimalRotation{N,T} # which is <: AbstractMatrix{T}
    mat::SMatrix{N, N, T, L} # The final parameter to SMatrix is the "length" of the matrix, 3 × 3 = 9
    InfinitesimalRotMatrix{N,T,L}(x::AbstractArray) where {N,T,L} = new{N,T,L}(convert(SMatrix{N,N,T,L}, x))
    # fixes #49 ambiguity introduced in StaticArrays 0.6.5
    InfinitesimalRotMatrix{N,T,L}(x::StaticArray) where {N,T,L} = new{N,T,L}(convert(SMatrix{N,N,T,L}, x))
end
InfinitesimalRotMatrix(x::SMatrix{N,N,T,L}) where {N,T,L} = InfinitesimalRotMatrix{N,T,L}(x)

Base.zero(::Type{InfinitesimalRotMatrix}) = error("The dimension of rotation is not specified.")

# These functions (plus size) are enough to satisfy the entire StaticArrays interface:
for N = 2:3
    L = N*N
    InfinitesimalRotMatrixN = Symbol(:InfinitesimalRotMatrix, N)
    @eval begin
        @inline InfinitesimalRotMatrix(t::NTuple{$L})  = InfinitesimalRotMatrix(SMatrix{$N,$N}(t))
        @inline (::Type{InfinitesimalRotMatrix{$N}})(t::NTuple{$L}) = InfinitesimalRotMatrix(SMatrix{$N,$N}(t))
        @inline InfinitesimalRotMatrix{$N,T}(t::NTuple{$L}) where {T} = InfinitesimalRotMatrix(SMatrix{$N,$N,T}(t))
        @inline InfinitesimalRotMatrix{$N,T,$L}(t::NTuple{$L}) where {T} = InfinitesimalRotMatrix(SMatrix{$N,$N,T}(t))
        const $InfinitesimalRotMatrixN{T} = InfinitesimalRotMatrix{$N, T, $L}
    end
end
Base.@propagate_inbounds Base.getindex(r::InfinitesimalRotMatrix, i::Int) = r.mat[i]
@inline Base.Tuple(r::InfinitesimalRotMatrix) = Tuple(r.mat)

@inline InfinitesimalRotMatrix(θ::Real) = InfinitesimalRotMatrix{2}(θ)
@inline function (::Type{InfinitesimalRotMatrix{2}})(θ::Real)
    InfinitesimalRotMatrix(@SMatrix [zero(θ) -θ; θ zero(θ)])
end
@inline function InfinitesimalRotMatrix{2,T}(θ::Real) where T
    InfinitesimalRotMatrix(@SMatrix T[zero(θ) -θ; θ zero(θ)])
end

# A rotation is more-or-less defined as being an orthogonal (or unitary) matrix
Base.:-(r::InfinitesimalRotMatrix) = InfinitesimalRotMatrix(-r.mat)

# By default, composition of rotations will go through InfinitesimalRotMatrix, unless overridden
@inline Base.:+(r1::InfinitesimalRotation, r2::InfinitesimalRotation) = InfinitesimalRotMatrix(r1) + InfinitesimalRotMatrix(r2)
@inline Base.:+(r1::InfinitesimalRotMatrix, r2::InfinitesimalRotation) = r1 + InfinitesimalRotMatrix(r2)
@inline Base.:+(r1::InfinitesimalRotation, r2::InfinitesimalRotMatrix) = InfinitesimalRotMatrix(r1) + r2
@inline Base.:+(r1::InfinitesimalRotMatrix, r2::InfinitesimalRotMatrix) = InfinitesimalRotMatrix(r1.mat + r2.mat)

@inline Base.:*(t::Number, r::InfinitesimalRotMatrix) = InfinitesimalRotMatrix(t*r.mat)
@inline Base.:*(r::InfinitesimalRotMatrix, t::Number) = t*r

# Special case multiplication of 2×2 rotation matrices: speedup using skew-symmetricity
@inline function Base.:+(r1::InfinitesimalRotMatrix{2}, r2::InfinitesimalRotMatrix{2})
    v = r1[2,1]+r2[2,1]
    s = @SMatrix [0  -v
                  v   0]
    return InfinitesimalRotMatrix(s)
end

# Special case multiplication of 3×3 rotation matrices: speedup using skew-symmetricity
@inline function Base.:+(r1::InfinitesimalRotMatrix{3}, r2::InfinitesimalRotMatrix{3})
    x = r1[6]+r2[6]
    y = r1[7]+r2[7]
    z = r1[2]+r2[2]
    s = @SMatrix [ 0 -z  y
                   z  0 -x
                  -y  x  0]

    return InfinitesimalRotMatrix(s)
end

"""
    struct InfinitesimalAngle2d{T} <: InfinitesimalRotation{2,T}
        v::T
    end

A 2×2 infinitesimal rotation matrix (i.e. skew-symmetric matrix).
[ 0 -v
  v  0 ]
"""
struct InfinitesimalAngle2d{T} <: InfinitesimalRotation{2,T}
    v::T
    InfinitesimalAngle2d{T}(r) where T = new{T}(r)
end

@inline function InfinitesimalAngle2d(r::T) where T <: Number
    InfinitesimalAngle2d{T}(r)
end

@inline InfinitesimalAngle2d(r::InfinitesimalRotation{2}) = InfinitesimalAngle2d(r[2,1])
@inline InfinitesimalAngle2d{T}(r::InfinitesimalRotation{2}) where {T} = InfinitesimalAngle2d{T}(r[2,1])

@inline Base.zero(::Type{A}) where {A<: InfinitesimalAngle2d} = A(0.0)

@inline Base.:+(r1::InfinitesimalAngle2d, r2::InfinitesimalAngle2d) = InfinitesimalAngle2d(r1.v + r2.v)
@inline Base.:-(r1::InfinitesimalAngle2d, r2::InfinitesimalAngle2d) = InfinitesimalAngle2d(r1.v - r2.v)
@inline Base.:*(t::Number, r::InfinitesimalAngle2d) = InfinitesimalAngle2d(t*r.v)
@inline Base.:*(r::InfinitesimalAngle2d, t::Number) = t*r
@inline Base.:-(r::InfinitesimalAngle2d) = InfinitesimalAngle2d(-r.v)

@inline function Base.getindex(r::InfinitesimalAngle2d{T}, i::Int) where T
    if i == 1
        zero(T)
    elseif i == 2
        r.v
    elseif i == 3
        -r.v
    elseif i == 4
        zero(T)
    else
        throw(BoundsError(r,i))
    end
end


"""
    struct InfinitesimalRotationVec{T} <: InfinitesimalRotation{2,T}
        x::T
        y::T
        z::T
    end

A 3×3 infinitesimal rotation matrix (i.e. skew-symmetric matrix).
[ 0 -z  y
  z  0 -x
 -y  x  0]
"""
struct InfinitesimalRotationVec{T} <: InfinitesimalRotation{3,T}
    x::T
    y::T
    z::T
    InfinitesimalRotationVec{T}(x,y,z) where T = new{T}(x,y,z)
end

@inline function InfinitesimalRotationVec(x::X,y::Y,z::Z) where {X<:Number, Y<:Number, Z<:Number}
    InfinitesimalRotationVec{promote_type(X,Y,Z)}(x, y, z)
end

@inline InfinitesimalRotationVec(r::InfinitesimalRotation{3}) = InfinitesimalRotationVec(r[6], r[7], r[2])
@inline InfinitesimalRotationVec{T}(r::InfinitesimalRotation{3}) where {T} = InfinitesimalRotationVec{T}(r[6], r[7], r[2])

@inline Base.zero(::Type{R}) where {R<: InfinitesimalRotationVec} = R(0.0,0.0,0.0)

@inline Base.:+(r1::InfinitesimalRotationVec, r2::InfinitesimalRotationVec) = InfinitesimalRotationVec(r1.x + r2.x, r1.y + r2.y, r1.z + r2.z)
@inline Base.:-(r1::InfinitesimalRotationVec, r2::InfinitesimalRotationVec) = InfinitesimalRotationVec(r1.x - r2.x, r1.y - r2.y, r1.z - r2.z)
@inline Base.:*(t::Number, r::InfinitesimalRotationVec) = InfinitesimalRotationVec(t*r.x, t*r.y, t*r.z)
@inline Base.:*(r::InfinitesimalRotationVec, t::Number) = t*r
@inline Base.:-(r::InfinitesimalRotationVec) = InfinitesimalRotationVec(-r.x, -r.y, -r.z)

@inline function Base.getindex(r::InfinitesimalRotationVec{T}, i::Int) where T
    if i == 1
        zero(T)
    elseif i == 2
        r.z
    elseif i == 3
        -r.y
    elseif i == 4
        -r.z
    elseif i == 5
        zero(T)
    elseif i == 6
        r.x
    elseif i == 7
        r.y
    elseif i == 8
        -r.x
    elseif i == 9
        zero(T)
    else
        throw(BoundsError(r,i))
    end
end


################################################################################
################################################################################


# A simplification and specialization of the Base.show function for AbstractArray makes
# everything sensible at the REPL.
function Base.show(io::IO, ::MIME"text/plain", X::InfinitesimalRotation)
    if !haskey(io, :compact)
        io = IOContext(io, :compact => true)
    end
    summary(io, X)
    if !isa(X, InfinitesimalRotMatrix)
        n_fields = length(fieldnames(typeof(X)))
        print(io, "(")
        for i = 1:n_fields
            print(io, getfield(X, i))
            if i < n_fields
                print(io, ", ")
            end
        end
        print(io, ")")
    end
    print(io, ":")
    println(io)
    io = IOContext(io, :typeinfo => eltype(X))
    Base.print_array(io, X)
end
