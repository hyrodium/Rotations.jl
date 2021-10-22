# Quaternion and Related Parameters

## `UnitQuaternion`
A 3D rotation parameterized by a unit quaternion.
Note that the constructor will renormalize the quaternion to be a unit quaternion, and that although they follow the same multiplicative *algebra* as quaternions, it is better to think of `UnitQuaternion` as a 3×3 matrix rather than as a quaternion *number*.

## `RodriguesParams`
A 3-parameter representation of 3D rotations that has a singularity at 180 degrees.
They can be interpreted as a projection of the unit quaternion onto the plane tangent to the quaternion identity.
They are computationally efficient and do not have a sign ambiguity.

## `MRP` (Modified Rodrigues Parameters)

A 3D rotation encoded by the stereographic projection of a unit quaternion.
This projection can be visualized as a pin hole camera, with the pin hole matching the quaternion `[-1,0,0,0]` and the image plane containing the origin and having normal direction `[1,0,0,0]`.
The "null rotation" `Quaternion(1.0,0,0,0)` then maps to the `MRP(0,0,0)`

These are similar to the Rodrigues vector in that the axis direction is stored in an unnormalized form, and the rotation angle is encoded in the length of the axis.
This type has the nice property that the derivatives of the rotation matrix w.r.t. the `MRP` parameters are rational functions, making the `MRP` type a good choice for differentiation / optimization.

They are frequently used in aerospace applications since they are a 3-parameter representation whose singularity happens at 360 degrees.
In practice, the singularity can be avoided with some switching logic between one of two equivalent MRPs (obtained by projecting the negated quaternion).
