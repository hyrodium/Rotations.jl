# Rotations

3D rotations made easy in Julia

## Installation
```
pkg> add Rotations
```

## First example

```@repl
using Rotations, StaticArrays

# create the null rotation (identity matrix)
id = one(RotMatrix{3, Float64})

# create a random rotation matrix (uniformly distributed over all 3D rotations)
r = rand(RotMatrix{3}) # uses Float64 by default

# create a point
p = SVector(1.0, 2.0, 3.0) # from StaticArrays.jl, but could use any AbstractVector...

# convert to a quaternion (Quat) and rotate the point
q = UnitQuaternion(r)
p_rotated = q * p

# Compose rotations
q2 = rand(UnitQuaternion)
q3 = q * q2

# Take the inverse (equivalent to transpose)
q_inv = transpose(q)
q_inv == inv(q)
p ≈ q_inv * (q * p)
q4 = q3 / q2  # q4 = q3 * inv(q2)
q5 = q3 \ q2  # q5 = inv(q3) * q2

# convert to a Modified Rodrigues Parameter (aka Stereographic quaternion projection, recommended for applications with differentiation)
spq = MRP(r)

# convert to the Angle-axis parameterization, or related Rotation vector
aa = AngleAxis(r)
rv = RotationVec(r)
ϕ = rotation_angle(r)
v = rotation_axis(r)

# convert to Euler angles, composed of X/Y/Z axis rotations (Z applied first)
# (all combinations of "RotABC" are defined)
r_xyz = RotXYZ(r)

# Rotation about the X axis by 0.1 radians
r_x = RotX(0.1)

# Composing axis rotations together automatically results in Euler parameterization
RotX(0.1) * RotY(0.2) * RotZ(0.3) === RotXYZ(0.1, 0.2, 0.3)

# Can calculate Jacobian - derivatives of rotations with respect to parameters
j1 = Rotations.jacobian(RotMatrix, q) # How does the matrix change w.r.t the 4 Quat parameters?
j2 = Rotations.jacobian(q, p) # How does the rotated point q*p change w.r.t. the 4 Quat parameters?
# ... all Jacobian's involving RotMatrix, MRP and Quat are implemented
# (MRP is ideal for optimization purposes - no constaints/singularities)
```