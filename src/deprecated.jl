# Deprecate UnitQuaternion => QuatRotation
Base.@deprecate_binding UnitQuaternion QuatRotation true

Base.@deprecate rotation_between(from::AbstractVector, to::AbstractVector) rotation_between(SVector{3}(from), SVector{3}(to))
