## log
# 2d
Base.log(R::Angle2d) = Angle2dGenerator(R.theta)
Base.log(R::RotMatrix{2}) = RotMatrixGenerator(log(Angle2d(R)))
#= We can define log for Rotation{2} like this,
but the subtypes of Rotation{2} are only Angle2d and RotMatrix{2},
so we don't need this defnition. =#
# Base.log(R::Rotation{2}) = log(Angle2d(R))

# 3d
Base.log(R::RotationVec) = RotationVecGenerator(R.sx,R.sy,R.sz)
Base.log(R::RotMatrix{3}) = RotMatrixGenerator(log(RotationVec(R)))
Base.log(R::Rotation{3}) = log(RotationVec(R))

# General dimensions
function Base.log(R::RotMatrix{N}) where N
    # This will be faster when log(::SMatrix) is implemented in StaticArrays.jl
    @static if VERSION < v"1.7"
        # This if block is related to this PR.
        # https://github.com/JuliaLang/julia/pull/40573
        S = SMatrix{N,N}(real(log(Matrix(R))))
    else
        S = SMatrix{N,N}(log(Matrix(R)))
    end
    RotMatrixGenerator((S-S')/2)
end

## exp
# 2d
Base.exp(R::Angle2dGenerator) = Angle2d(R.v)
Base.exp(R::RotMatrixGenerator{2}) = RotMatrix(exp(Angle2dGenerator(R)))
# Same as log(R::Rotation{2}), this definition is not necessary until someone add a new subtype of RotationGenerator{2}.
# Base.exp(R::RotationGenerator{2}) = exp(Angle2dGenerator(R))

# 3d
Base.exp(R::RotationVecGenerator) = RotationVec(R.x,R.y,R.z)
Base.exp(R::RotMatrixGenerator{3}) = RotMatrix(exp(RotationVecGenerator(R)))
# Same as log(R::Rotation{2}), this definition is not necessary until someone add a new subtype of RotationGenerator{3}.
# Base.exp(R::RotationGenerator{3}) = exp(RotationVecGenerator(R))

# General dimensions
Base.exp(R::RotMatrixGenerator{N}) where N = RotMatrix(exp(SMatrix(R)))

## sqrt
# 2d
Base.sqrt(r::Angle2d) = Angle2d(r.theta/2)
Base.sqrt(r::RotMatrix{2}) = RotMatrix(sqrt(Angle2d(r)))

# 3d
Base.sqrt(r::RotX) = RotX(r.theta/2)
Base.sqrt(r::RotY) = RotY(r.theta/2)
Base.sqrt(r::RotZ) = RotZ(r.theta/2)
Base.sqrt(r::AngleAxis) = AngleAxis(r.theta/2, r.axis_x, r.axis_y, r.axis_z)
Base.sqrt(r::RotationVec) = RotaitonVec(r.sx/2, r.sy/2, r.sz/2)
Base.sqrt(r::QuatRotation) = QuatRotation(sqrt(r.q))
Base.sqrt(r::RotMatrix{3}) = RotMatrix{3}(sqrt(QuatRotation(r)))
Base.sqrt(r::RodriguesParam) = RodriguesParam(sqrt(QuatRotation(r)))
Base.sqrt(r::MRP) = MRP(sqrt(QuatRotation(r)))
Base.sqrt(r::Rotation{3}) = sqrt(QuatRotation(r))

# General dimensions
Base.sqrt(r::Rotation{N}) where N = RotMatrix(sqrt(r))

## cbrt
# 2d
Base.cbrt(r::Angle2d) = Angle2d(r.theta/3)
Base.cbrt(r::RotMatrix{2}) = RotMatrix(cbrt(Angle2d(r)))

# 3d
Base.cbrt(r::RotX) = RotX(r.theta/3)
Base.cbrt(r::RotY) = RotY(r.theta/3)
Base.cbrt(r::RotZ) = RotZ(r.theta/3)
Base.cbrt(r::AngleAxis) = AngleAxis(r.theta/3, r.axis_x, r.axis_y, r.axis_z)
Base.cbrt(r::RotationVec) = RotaitonVec(r.sx/3, r.sy/3, r.sz/3)
Base.cbrt(r::QuatRotation) = QuatRotation(cbrt(r.q))
Base.cbrt(r::RotMatrix{3}) = RotMatrix{3}(cbrt(QuatRotation(r)))
Base.cbrt(r::RodriguesParam) = RodriguesParam(cbrt(QuatRotation(r)))
Base.cbrt(r::MRP) = MRP(cbrt(QuatRotation(r)))
Base.cbrt(r::Rotation{3}) = cbrt(QuatRotation(r))

# General dimensions
Base.cbrt(r::Rotation{N}) where N = exp(log(r)/3)
