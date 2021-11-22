## log
# 3d
function Base.log(R::RotationVec)
    x, y, z = params(R)
    return RotationVecGenerator(x,y,z)
end

function Base.log(R::Rotation{3})
    log(RotationVec(R))
end


# 2d
function Base.log(R::Angle2d)
    θ, = params(R)
    return Angle2dGenerator(θ)
end

function Base.log(R::Rotation{2})
    log(Angle2d(R))
end


## exp
# 3d
function Base.exp(R::RotationVecGenerator)
    return RotationVec(R.x,R.y,R.z)
end

function Base.exp(R::RotationGenerator{3})
    exp(RotationVecGenerator(R))
end


# 2d
function Base.exp(R::Angle2dGenerator)
    return Angle2d(R.v)
end

function Base.exp(R::RotationGenerator{2})
    exp(Angle2dGenerator(R))
end
