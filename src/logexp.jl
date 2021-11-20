## log
# 3d
function Base.log(R::RotationVec)
    x, y, z = params(R)
    return InfinitesimalRotationVec(x,y,z)
end

function Base.log(R::Rotation{3})
    log(RotationVec(R))
end


# 2d
function Base.log(R::Angle2d)
    θ, = params(R)
    return InfinitesimalAngle2d(θ)
end

function Base.log(R::Rotation{2})
    log(Angle2d(R))
end


## exp
# 3d
function Base.exp(R::InfinitesimalRotationVec)
    return RotationVec(R.x,R.y,R.z)
end

function Base.exp(R::InfinitesimalRotation{3})
    exp(InfinitesimalRotationVec(R))
end


# 2d
function Base.exp(R::InfinitesimalAngle2d)
    return Angle2d(R.v)
end

function Base.exp(R::InfinitesimalRotation{2})
    exp(InfinitesimalAngle2d(R))
end
