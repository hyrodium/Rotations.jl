# 2D Rotation Matrix

An ``2 \times 2`` rotation matrix storing the rotation.
This is a simple wrapper for a [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) `SMatrix{2,2,T}`.
A rotation matrix ``R`` should have the property ``I = R R^\top``, but this isn't enforced by the constructor.
On the other hand, all the types below are guaranteed to be "proper" rotations for all input parameters (equivalently: parity conserving, in ``SO(2)``, ``\det(R) = 1``, or a rotation without reflection).

```@setup rotmatrix2
using Rotations
```

## `RotMatrix2`

**example**

```@repl rotmatrix2
t = 1.2  # rotation angle
m = [cos(t) -sin(t);sin(t) cos(t)]
RotMatrix{2}(m)  # construct from matrix
RotMatrix{2}(t)  # construct from angle
```
