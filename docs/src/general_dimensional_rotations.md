# General dimensional Rotation

Our main goal is to efficiently represent rotations in lower dimensions, such as 2D and 3D, but some operations on rotations in general dimensions are also supported.

```@setup rotmatrixN
using Rotations
using StaticArrays
```

**example**

```@repl rotmatrixN
r1 = one(RotMatrix{4})  # generate identity rotation matrix
m = @SMatrix rand(4,4)
r2 = nearest_rotation(m)  # nearest rotation matrix from given matrix
r3 = rand(RotMatrix{4})  # random rotation in SO(4)
r1*r2/r3  # multiplication and division
s = log(r2)  # logarithm of RotMatrix is a RotMatrixGenerator
exp(s)  # exponential of RotMatrixGenerator is a RotMatrix
```
