
# Deprecate RodriguesVec => RotationVec
Base.@deprecate_binding RodriguesVec RotationVec true

# Deprecate Quat => QuatRotation
Base.@deprecate_binding Quat QuatRotation true

# Deprecate SPQuat => MRP
Base.@deprecate_binding SPQuat MRP true
