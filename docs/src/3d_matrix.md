# 3D Rotation Matrix

An ``3 \times 3`` rotation matrix storing the rotation.
This is a simple wrapper for a [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) `SMatrix{3,3,T}`.
A rotation matrix ``R`` should have the property ``I = R R^\top``, but this isn't enforced by the constructor.
On the other hand, all the types below are guaranteed to be "proper" rotations for all input parameters (equivalently: parity conserving, in ``SO(3)``, ``\det(R) = 1``, or a rotation without reflection).

```@setup rotmatrix3
using Rotations
```

## `RotMatrix3`

**example**

```@repl rotmatrix3
t = 1.2
m1 = [cos(t) -sin(t) 0;sin(t) cos(t) 0;0 0 1]
RotMatrix{3}(m1)  # construct from matrix
m2 = RotZ(t)  # Other rotation type
RotMatrix{3}(m2)  # convert m2 to RotMatrix
RotMatrix(m2)  # Type parameter can be omitted.
```
