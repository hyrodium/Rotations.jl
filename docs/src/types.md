# Types in Rotations.jl

## Type hierarchy

```@setup abc
using InteractiveUtils
```

```@repl abc
using Rotations, StaticArrays
Rotation <: StaticMatrix <: AbstractMatrix
subtypes(Rotation{2})
subtypes(Rotation{3})
```

## Abstract rotations
If a matrix ``R`` satisfies ``R^\top = R^{-1}`` and ``\det(R)=1``, the matrix ``R`` is said to be *rotation matrix*.

In `Rotations.jl`, there's an abstract type for rotations matrix, `Rotation{D}`.
Where `D` is a size of the rotation matrix.

## 2D rotations
* `RotMatrix{2, T, L} where {T, L}`
    * Rotation matrix in 2 dim.
* `Angle2d`
    * Parametrized with rotational angle.

## 3D rotations
* `RotMatrix{3, T, L} where {T, L}`
    * Rotation matrix in 3 dim.
* `RotX`, `RotYZ`, `RotXYZ` and etc.
    * Euler angles.
* `AngleAxis`
    * Rotation around given axis and angle.
* `RotationVec`
    * Rotation around given axis. The length of axis vector represents its angle.
* `UnitQuaternion`
    * A 3D rotation parameterized by a unit quaternion.
* `MRP`
    * A 3D rotation encoded by the stereographic projection of a unit quaternion.
