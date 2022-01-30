## log
# 2d
Base.log(R::Angle2d) = Angle2dGenerator(R.theta)
Base.log(R::Rotation{2}) = log(Angle2d(R))
Base.log(R::RotMatrix{2}) = RotMatrixGenerator(log(Angle2d(R)))

# 3d
Base.log(R::RotationVec) = RotationVecGenerator(R.sx,R.sy,R.sz)
Base.log(R::Rotation{3}) = log(RotationVec(R))
Base.log(R::RotMatrix{3}) = RotMatrixGenerator(log(RotationVec(R)))

# General dimensions
# This will be faster when log(::SMatrix) is implemented in StaticArrays.jl
Base.log(R::RotMatrix{N}) where N = RotMatrixGenerator(SMatrix{N,N}(log(Matrix(R))))

## exp
# 2d
Base.exp(R::Angle2dGenerator) = Angle2d(R.v)
Base.exp(R::RotationGenerator{2}) = exp(Angle2dGenerator(R))
Base.exp(R::RotMatrixGenerator{2}) = RotMatrix(exp(Angle2dGenerator(R)))

# 3d
Base.exp(R::RotationVecGenerator) = RotationVec(R.x,R.y,R.z)
Base.exp(R::RotationGenerator{3}) = exp(RotationVecGenerator(R))
Base.exp(R::RotMatrixGenerator{3}) = RotMatrix(exp(RotationVecGenerator(R)))

# General dimensions
Base.exp(R::RotMatrixGenerator{N}) where N = RotMatrix(exp(SMatrix(R)))
